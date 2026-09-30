//! Per-cell metadata (sample, condition, cell type, cluster, ...) used for
//! stratified references and provenance.
//!
//! Sources:
//! - 10x directories: `metadata.tsv[.gz]` next to the matrix (or `--metadata`):
//!   a header line, first column = barcode, tab-separated.
//! - AnnData: `obs/<column>` string datasets or categorical groups
//!   (`codes` + `categories`, anndata >= 0.8) and the legacy
//!   `obs/__categories/<column>` layout.
//!
//! Columns are re-ordered to the canonical (sorted) cell order of the main
//! matrix; cells missing from the metadata get an empty string.

use std::collections::{BTreeMap, HashMap};
use std::io::BufRead;
use std::path::Path;

use tracing::{info, warn};

use crate::input::error::InputError;

/// Column aliases, in priority order, for the roles the pipeline understands.
pub const CELL_TYPE_ALIASES: &[&str] = &[
    "cell_type",
    "celltype",
    "cell_type_annotation",
    "annotation",
    "cell_ontology_class",
    "predicted.celltype",
];
pub const CLUSTER_ALIASES: &[&str] =
    &["cluster", "clusters", "leiden", "louvain", "seurat_clusters", "cluster_id"];
pub const SAMPLE_ALIASES: &[&str] = &["sample", "sample_id", "orig.ident", "batch", "library"];
pub const CONDITION_ALIASES: &[&str] = &["condition", "group", "treatment", "disease"];
/// Boolean-like doublet calls from upstream tools (Scrublet, scDblFinder,
/// DoubletFinder); truthy values: `true`, `1`, `yes`, `doublet`.
pub const DOUBLET_ALIASES: &[&str] = &[
    "predicted_doublet",
    "doublet",
    "is_doublet",
    "scDblFinder.class",
    "scdblfinder_class",
    "doublet_class",
    "DF.classifications",
];

/// Truthy interpretation of a doublet-call cell value.
pub fn is_doublet_value(value: &str) -> bool {
    matches!(
        value.trim().to_ascii_lowercase().as_str(),
        "true" | "1" | "yes" | "doublet"
    )
}

#[derive(Debug, Clone, Default)]
pub struct CellMetadata {
    /// Column name -> one value per canonical cell (empty when unknown).
    pub columns: BTreeMap<String, Vec<String>>,
    pub source: String,
    /// Canonical cells that had no metadata row.
    pub cells_without_metadata: usize,
}

impl CellMetadata {
    pub fn column(&self, name: &str) -> Option<&[String]> {
        self.columns.get(name).map(|v| v.as_slice())
    }

    /// First alias that exists as a column (case-insensitive), with its values.
    pub fn resolve(&self, aliases: &[&str]) -> Option<(&str, &[String])> {
        for alias in aliases {
            if let Some((name, values)) = self
                .columns
                .iter()
                .find(|(name, _)| name.eq_ignore_ascii_case(alias))
            {
                return Some((name.as_str(), values.as_slice()));
            }
        }
        None
    }

    pub fn is_empty(&self) -> bool {
        self.columns.is_empty()
    }
}

/// Reads a `metadata.tsv[.gz]` table and aligns it to `cell_names`.
pub fn read_metadata_tsv(path: &Path, cell_names: &[String]) -> Result<CellMetadata, InputError> {
    let reader = kira_scio::open_maybe_gz_existing(path)
        .map_err(|e| InputError::UnsupportedInput(e.message))?;
    let mut lines = reader.lines();
    let header = match lines.next() {
        Some(h) => h.map_err(|e| InputError::io(path, e))?,
        None => return Err(InputError::InvalidMetadata(format!("{}: empty file", path.display()))),
    };
    let header: Vec<String> = header
        .trim_end_matches(['\r', '\n'])
        .split('\t')
        .map(|s| s.trim().to_string())
        .collect();
    if header.len() < 2 {
        return Err(InputError::InvalidMetadata(format!(
            "{}: expected a barcode column and at least one metadata column",
            path.display()
        )));
    }

    let mut rows: HashMap<String, Vec<String>> = HashMap::new();
    for (line_no, line) in lines.enumerate() {
        let line = line.map_err(|e| InputError::io(path, e))?;
        let t = line.trim_end_matches(['\r', '\n']);
        if t.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = t.split('\t').collect();
        if fields.len() != header.len() {
            return Err(InputError::InvalidMetadata(format!(
                "{}:{}: {} fields, header has {}",
                path.display(),
                line_no + 2,
                fields.len(),
                header.len()
            )));
        }
        rows.insert(
            fields[0].trim().to_string(),
            fields[1..].iter().map(|f| f.trim().to_string()).collect(),
        );
    }

    let mut columns: BTreeMap<String, Vec<String>> = header[1..]
        .iter()
        .map(|name| (name.clone(), Vec::with_capacity(cell_names.len())))
        .collect();
    let mut missing = 0usize;
    for cell in cell_names {
        match rows.get(cell) {
            Some(values) => {
                for (name, value) in header[1..].iter().zip(values) {
                    columns.get_mut(name).unwrap().push(value.clone());
                }
            }
            None => {
                missing += 1;
                for name in &header[1..] {
                    columns.get_mut(name).unwrap().push(String::new());
                }
            }
        }
    }
    if missing > 0 {
        warn!(
            path = %path.display(),
            cells_without_metadata = missing,
            "cells without a metadata row"
        );
    }
    info!(path = %path.display(), columns = columns.len(), "cell metadata loaded");
    Ok(CellMetadata {
        columns,
        source: format!("tsv:{}", path.display()),
        cells_without_metadata: missing,
    })
}

/// Reads string / categorical `obs` columns of an AnnData file.
/// `obs_order` is the file's own cell order (obs index), `cell_names` the
/// canonical order to align to.
pub fn read_metadata_h5ad(
    path: &Path,
    obs_order: &[String],
    cell_names: &[String],
) -> Result<CellMetadata, InputError> {
    let file = hdf5::File::open(path).map_err(|e| InputError::UnsupportedInput(e.to_string()))?;
    let obs = match file.group("obs") {
        Ok(g) => g,
        Err(_) => return Ok(CellMetadata::default()),
    };

    let mut raw: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for name in obs.member_names().unwrap_or_default() {
        if name.starts_with('_') || name == "__categories" {
            continue;
        }
        if let Some(values) = read_obs_column(&obs, &name)
            && values.len() == obs_order.len()
        {
            raw.insert(name, values);
        }
    }
    if raw.is_empty() {
        return Ok(CellMetadata::default());
    }

    let position: HashMap<&str, usize> = obs_order
        .iter()
        .enumerate()
        .map(|(i, b)| (b.as_str(), i))
        .collect();
    let mut missing = 0usize;
    let mut columns = BTreeMap::new();
    for (name, values) in raw {
        let aligned: Vec<String> = cell_names
            .iter()
            .map(|c| position.get(c.as_str()).map_or(String::new(), |&i| values[i].clone()))
            .collect();
        columns.insert(name, aligned);
    }
    for c in cell_names {
        if !position.contains_key(c.as_str()) {
            missing += 1;
        }
    }
    info!(path = %path.display(), columns = columns.len(), "cell metadata loaded from obs");
    Ok(CellMetadata {
        columns,
        source: format!("h5ad-obs:{}", path.display()),
        cells_without_metadata: missing,
    })
}

fn read_obs_column(obs: &hdf5::Group, name: &str) -> Option<Vec<String>> {
    // anndata >= 0.8 categorical: group with `codes` and `categories`.
    if let Ok(group) = obs.group(name) {
        let categories = read_strings(&group.dataset("categories").ok()?)?;
        let codes = read_codes(&group.dataset("codes").ok()?)?;
        return Some(decode(&codes, &categories));
    }
    let ds = obs.dataset(name).ok()?;
    if let Some(values) = read_strings(&ds) {
        return Some(values);
    }
    // Legacy categorical: integer codes + obs/__categories/<name>.
    if let Ok(cats) = obs.dataset(&format!("__categories/{name}"))
        && let Some(categories) = read_strings(&cats)
        && let Some(codes) = read_codes(&ds)
    {
        return Some(decode(&codes, &categories));
    }
    None
}

fn decode(codes: &[i64], categories: &[String]) -> Vec<String> {
    codes
        .iter()
        .map(|&c| {
            if c < 0 {
                String::new()
            } else {
                categories.get(c as usize).cloned().unwrap_or_default()
            }
        })
        .collect()
}

fn read_codes(ds: &hdf5::Dataset) -> Option<Vec<i64>> {
    ds.read_raw::<i64>()
        .or_else(|_| ds.read_raw::<i32>().map(|v| v.into_iter().map(i64::from).collect()))
        .or_else(|_| ds.read_raw::<i16>().map(|v| v.into_iter().map(i64::from).collect()))
        .or_else(|_| ds.read_raw::<i8>().map(|v| v.into_iter().map(i64::from).collect()))
        .ok()
}

fn read_strings(ds: &hdf5::Dataset) -> Option<Vec<String>> {
    use hdf5::types::{VarLenAscii, VarLenUnicode};
    if let Ok(v) = ds.read_raw::<VarLenUnicode>() {
        return Some(v.into_iter().map(|s| s.to_string()).collect());
    }
    if let Ok(v) = ds.read_raw::<VarLenAscii>() {
        return Some(v.into_iter().map(|s| s.to_string()).collect());
    }
    None
}
