//! Readers for spliced / unspliced / ambiguous count layers (input level L1).
//!
//! Supported sources:
//! - a directory with `spliced.mtx[.gz]`, `unspliced.mtx[.gz]` and optionally
//!   `ambiguous.mtx[.gz]` (STARsolo `Velocyto/{raw,filtered}`, kb-python
//!   `--workflow lamanno`/`nac` outputs). When the directory carries its own
//!   `barcodes.tsv[.gz]`, cells are matched to the main matrix by barcode;
//!   otherwise the layers must have the main matrix's dimensions and order.
//! - AnnData `layers/spliced`, `layers/unspliced` (and `layers/ambiguous`)
//!   sparse groups, which share `X`'s row/column order.
//!
//! Layers are returned as triplets in the main matrix's *original* gene and
//! cell index space; stage 1 reindexes them together with the main matrix.

use std::collections::HashMap;
use std::io::BufRead;
use std::path::{Path, PathBuf};

use tracing::{info, warn};

use crate::input::error::InputError;

/// Where the layers of a dataset live.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LayerLocation {
    /// Directory with `spliced.mtx`/`unspliced.mtx` (+ optional barcodes/features).
    MtxDir(PathBuf),
    /// `layers/` groups inside an AnnData file.
    H5ad(PathBuf),
}

impl LayerLocation {
    pub fn describe(&self) -> String {
        match self {
            LayerLocation::MtxDir(p) => format!("mtx-dir:{}", p.display()),
            LayerLocation::H5ad(p) => format!("h5ad-layers:{}", p.display()),
        }
    }
}

/// One layer as `(gene, cell, count)` triplets in the main matrix's index space.
#[derive(Debug, Clone, Default)]
pub struct RawLayer {
    pub triplets: Vec<(u32, u32, u32)>,
}

#[derive(Debug, Clone)]
pub struct RawLayers {
    pub spliced: RawLayer,
    pub unspliced: RawLayer,
    pub ambiguous: Option<RawLayer>,
    pub source: String,
    /// Main-matrix cells with no counterpart in the layer source.
    pub cells_without_layers: usize,
}

const SPLICED: &str = "spliced.mtx";
const UNSPLICED: &str = "unspliced.mtx";
const AMBIGUOUS: &str = "ambiguous.mtx";

/// Finds layer files for a 10x-style input directory.
///
/// Search order: an explicit override, the input directory itself, then the
/// STARsolo sibling layout (`<root>/Gene/<subset>` -> `<root>/Velocyto/<subset>`).
pub fn detect_mtx_layers(input_dir: &Path, override_path: Option<&Path>) -> Option<LayerLocation> {
    if let Some(p) = override_path {
        return Some(if p.is_dir() {
            LayerLocation::MtxDir(p.to_path_buf())
        } else {
            LayerLocation::H5ad(p.to_path_buf())
        });
    }
    if has_layer_files(input_dir) {
        return Some(LayerLocation::MtxDir(input_dir.to_path_buf()));
    }
    // STARsolo: Solo.out/Gene/raw <-> Solo.out/Velocyto/raw
    if let (Some(subset), Some(gene_dir)) = (input_dir.file_name(), input_dir.parent())
        && gene_dir
            .file_name()
            .is_some_and(|n| n == "Gene" || n == "GeneFull")
        && let Some(root) = gene_dir.parent()
    {
        let candidate = root.join("Velocyto").join(subset);
        if has_layer_files(&candidate) {
            return Some(LayerLocation::MtxDir(candidate));
        }
    }
    None
}

fn has_layer_files(dir: &Path) -> bool {
    dir.is_dir()
        && kira_scio::exists_plain_or_gz(&dir.join(SPLICED))
        && kira_scio::exists_plain_or_gz(&dir.join(UNSPLICED))
}

fn resolve(dir: &Path, base: &str) -> Option<PathBuf> {
    let plain = dir.join(base);
    if plain.is_file() {
        return Some(plain);
    }
    let gz = kira_scio::gz_path(&plain);
    if gz.is_file() { Some(gz) } else { None }
}

/// Reads the layers of `location` for a main matrix with `genes` x `barcodes`
/// (in the main matrix's original order).
pub fn read_layers(
    location: &LayerLocation,
    genes: &[String],
    barcodes: &[String],
) -> Result<RawLayers, InputError> {
    match location {
        LayerLocation::MtxDir(dir) => read_mtx_dir(dir, genes, barcodes),
        LayerLocation::H5ad(path) => read_h5ad_layers(path, genes.len(), barcodes.len()),
    }
}

// ---------------------------------------------------------------------------
// MTX directory
// ---------------------------------------------------------------------------

fn read_mtx_dir(
    dir: &Path,
    genes: &[String],
    barcodes: &[String],
) -> Result<RawLayers, InputError> {
    let spliced_path = resolve(dir, SPLICED)
        .ok_or_else(|| InputError::MissingFile(dir.join(SPLICED).display().to_string()))?;
    let unspliced_path = resolve(dir, UNSPLICED)
        .ok_or_else(|| InputError::MissingFile(dir.join(UNSPLICED).display().to_string()))?;
    let ambiguous_path = resolve(dir, AMBIGUOUS);

    // Cell mapping: by barcode when the layer directory carries its own
    // barcodes, positional otherwise.
    let layer_barcodes = match resolve(dir, "barcodes.tsv") {
        Some(p) => Some(read_lines(&p)?),
        None => None,
    };
    let (cell_map, cells_without_layers) = match layer_barcodes {
        Some(ref lb) if lb != barcodes => {
            let index: HashMap<&str, u32> = barcodes
                .iter()
                .enumerate()
                .map(|(i, b)| (b.as_str(), i as u32))
                .collect();
            let map: Vec<Option<u32>> = lb.iter().map(|b| index.get(b.as_str()).copied()).collect();
            let matched = map.iter().filter(|m| m.is_some()).count();
            let missing = barcodes.len().saturating_sub(matched);
            if matched == 0 {
                return Err(InputError::LayerMismatch(format!(
                    "no barcode of {} matches the main matrix",
                    dir.display()
                )));
            }
            if missing > 0 {
                warn!(
                    layer_dir = %dir.display(),
                    cells_without_layers = missing,
                    "some main-matrix cells have no column in the layer files; their layer counts are 0"
                );
            }
            (Some(map), missing)
        }
        _ => (None, 0),
    };

    let n_layer_cells = cell_map.as_ref().map(|m| m.len()).unwrap_or(barcodes.len());
    let read = |path: &Path| -> Result<RawLayer, InputError> {
        let (rows, cols, triplets) = parse_mtx_counts(path)?;
        if rows != genes.len() {
            return Err(InputError::LayerMismatch(format!(
                "{}: {} genes, main matrix has {}",
                path.display(),
                rows,
                genes.len()
            )));
        }
        if cols != n_layer_cells {
            return Err(InputError::LayerMismatch(format!(
                "{}: {} cells, expected {}",
                path.display(),
                cols,
                n_layer_cells
            )));
        }
        let triplets = match &cell_map {
            None => triplets,
            Some(map) => triplets
                .into_iter()
                .filter_map(|(g, c, k)| map[c as usize].map(|mc| (g, mc, k)))
                .collect(),
        };
        Ok(RawLayer { triplets })
    };

    let spliced = read(&spliced_path)?;
    let unspliced = read(&unspliced_path)?;
    let ambiguous = match ambiguous_path {
        Some(p) => Some(read(&p)?),
        None => None,
    };
    info!(
        layer_dir = %dir.display(),
        spliced_nnz = spliced.triplets.len(),
        unspliced_nnz = unspliced.triplets.len(),
        ambiguous = ambiguous.is_some(),
        "spliced/unspliced layers loaded"
    );
    Ok(RawLayers {
        spliced,
        unspliced,
        ambiguous,
        source: format!("mtx-dir:{}", dir.display()),
        cells_without_layers,
    })
}

fn read_lines(path: &Path) -> Result<Vec<String>, InputError> {
    let reader = kira_scio::open_maybe_gz_existing(path)
        .map_err(|e| InputError::UnsupportedInput(e.message))?;
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| InputError::io(path, e))?;
        let t = line.trim();
        if !t.is_empty() {
            // Only the first column matters (barcode files are single-column).
            out.push(t.split('\t').next().unwrap_or(t).to_string());
        }
    }
    Ok(out)
}

/// Parses a coordinate MatrixMarket file into `(n_rows, n_cols, (row, col, count))`
/// with 0-based indices. Values must be non-negative integers (a tolerance of
/// 1e-4 absorbs float round-trips).
/// `(n_rows, n_cols, (row, col, count) triplets)` of a parsed MatrixMarket file.
type ParsedCounts = (usize, usize, Vec<(u32, u32, u32)>);

pub(crate) fn parse_mtx_counts(path: &Path) -> Result<ParsedCounts, InputError> {
    const FRAC_TOL: f64 = 1e-4;
    let reader = kira_scio::open_maybe_gz_existing(path)
        .map_err(|e| InputError::UnsupportedInput(e.message))?;
    let mut header: Option<(usize, usize)> = None;
    let mut triplets = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| InputError::io(path, e))?;
        let t = line.trim();
        if t.is_empty() || t.starts_with('%') {
            continue;
        }
        let mut parts = t.split_whitespace();
        if header.is_none() {
            let r = parts.next().and_then(|s| s.parse::<usize>().ok());
            let c = parts.next().and_then(|s| s.parse::<usize>().ok());
            match (r, c) {
                (Some(r), Some(c)) => header = Some((r, c)),
                _ => return Err(InputError::InvalidMatrixMarket),
            }
            continue;
        }
        let (rows, cols) = header.unwrap();
        let row = parts.next().and_then(|s| s.parse::<usize>().ok());
        let col = parts.next().and_then(|s| s.parse::<usize>().ok());
        let val = parts.next().and_then(|s| s.parse::<f64>().ok());
        let (Some(row), Some(col), Some(val)) = (row, col, val) else {
            return Err(InputError::InvalidMatrixMarket);
        };
        if row == 0 || col == 0 || row > rows || col > cols {
            return Err(InputError::InvalidSparseMatrix);
        }
        if !val.is_finite() || val < 0.0 || val > u32::MAX as f64 || val.fract().abs() > FRAC_TOL {
            return Err(InputError::InvalidSparseMatrix);
        }
        let k = val.round() as u32;
        if k > 0 {
            triplets.push(((row - 1) as u32, (col - 1) as u32, k));
        }
    }
    let (rows, cols) = header.ok_or(InputError::InvalidMatrixMarket)?;
    Ok((rows, cols, triplets))
}

// ---------------------------------------------------------------------------
// AnnData layers
// ---------------------------------------------------------------------------

/// True when the file carries both `layers/spliced` and `layers/unspliced`.
pub fn h5ad_has_layers(path: &Path) -> bool {
    hdf5::File::open(path)
        .map(|f| f.group("layers/spliced").is_ok() && f.group("layers/unspliced").is_ok())
        .unwrap_or(false)
}

fn read_h5ad_layers(path: &Path, n_genes: usize, n_cells: usize) -> Result<RawLayers, InputError> {
    let file = hdf5::File::open(path).map_err(|e| InputError::UnsupportedInput(e.to_string()))?;
    let spliced = read_h5ad_sparse_layer(&file, "layers/spliced", n_genes, n_cells)?;
    let unspliced = read_h5ad_sparse_layer(&file, "layers/unspliced", n_genes, n_cells)?;
    let ambiguous = if file.group("layers/ambiguous").is_ok() {
        Some(read_h5ad_sparse_layer(
            &file,
            "layers/ambiguous",
            n_genes,
            n_cells,
        )?)
    } else {
        None
    };
    info!(
        path = %path.display(),
        spliced_nnz = spliced.triplets.len(),
        unspliced_nnz = unspliced.triplets.len(),
        ambiguous = ambiguous.is_some(),
        "spliced/unspliced layers loaded"
    );
    Ok(RawLayers {
        spliced,
        unspliced,
        ambiguous,
        source: format!("h5ad-layers:{}", path.display()),
        cells_without_layers: 0,
    })
}

fn read_h5ad_sparse_layer(
    file: &hdf5::File,
    group_path: &str,
    n_genes: usize,
    n_cells: usize,
) -> Result<RawLayer, InputError> {
    const FRAC_TOL: f32 = 1e-4;
    let group = file
        .group(group_path)
        .map_err(|_| InputError::MissingDataset(group_path.to_string()))?;
    let encoding =
        read_attr_string(&group, "encoding-type").unwrap_or_else(|| "csr_matrix".to_string());
    let indptr = read_i64_dataset(&group, "indptr")?;
    let indices = read_i64_dataset(&group, "indices")?;
    let data: Vec<f32> = group
        .dataset("data")
        .and_then(|d| d.read_raw::<f32>())
        .map_err(|e| InputError::UnsupportedInput(format!("{group_path}/data: {e}")))?;
    if indices.len() != data.len() {
        return Err(InputError::LayerMismatch(format!(
            "{group_path}: indices/data length mismatch"
        )));
    }

    // AnnData: CSR = rows are cells; CSC = columns are cells.
    let (n_major, major_is_cell) = match encoding.as_str() {
        "csr_matrix" => (n_cells, true),
        "csc_matrix" => (n_genes, false),
        other => {
            return Err(InputError::UnsupportedInput(format!(
                "{group_path}: unsupported encoding {other}"
            )));
        }
    };
    if indptr.len() != n_major + 1 {
        return Err(InputError::LayerMismatch(format!(
            "{group_path}: indptr length {} does not match the main matrix",
            indptr.len()
        )));
    }

    let mut triplets = Vec::with_capacity(data.len());
    for major in 0..n_major {
        let s = indptr[major] as usize;
        let e = indptr[major + 1] as usize;
        for idx in s..e {
            let minor = indices[idx] as usize;
            let v = data[idx];
            if !v.is_finite() || v < 0.0 || v.fract().abs() > FRAC_TOL {
                return Err(InputError::InvalidSparseMatrix);
            }
            let (gene, cell) = if major_is_cell {
                (minor, major)
            } else {
                (major, minor)
            };
            if gene >= n_genes || cell >= n_cells {
                return Err(InputError::LayerMismatch(format!(
                    "{group_path}: index out of range for the main matrix"
                )));
            }
            let k = v.round() as u32;
            if k > 0 {
                triplets.push((gene as u32, cell as u32, k));
            }
        }
    }
    Ok(RawLayer { triplets })
}

fn read_i64_dataset(group: &hdf5::Group, name: &str) -> Result<Vec<i64>, InputError> {
    let ds = group
        .dataset(name)
        .map_err(|_| InputError::MissingDataset(name.to_string()))?;
    ds.read_raw::<i64>()
        .or_else(|_| {
            ds.read_raw::<i32>()
                .map(|v| v.into_iter().map(i64::from).collect())
        })
        .map_err(|e| InputError::UnsupportedInput(format!("{name}: {e}")))
}

fn read_attr_string(group: &hdf5::Group, name: &str) -> Option<String> {
    use hdf5::types::VarLenUnicode;
    let attr = group.attr(name).ok()?;
    let value: VarLenUnicode = attr.read_scalar().ok()?;
    Some(value.to_string())
}
