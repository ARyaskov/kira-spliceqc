use std::path::Path;

use tracing::info;

use crate::expression::cache_writer::{CacheData, write_expr_bin};
use crate::expression::index::build_index;
use crate::expression::ExpressionMatrix;
use crate::expression::layers::{LayerMatrix, SplicedUnspliced};
use crate::expression::mmap::MmapExpressionMatrix;
use crate::expression::junctions::JunctionSet;
use crate::io::junctions::{RawJunctions, read_junctions};
use crate::io::layers::{RawLayer, RawLayers, read_layers};
use crate::input::error::InputError;
use crate::input::metadata::{CellMetadata, read_metadata_h5ad, read_metadata_tsv};
use crate::input::{InputDescriptor, InputKind};
use crate::io::{h5ad, mtx};

/// Stage 1 result: the mmap'd main matrix plus, when the input carries them,
/// the spliced/unspliced layers reindexed to the same gene/cell order.
#[derive(Debug)]
pub struct Stage1Output {
    pub matrix: MmapExpressionMatrix,
    pub layers: Option<SplicedUnspliced>,
    /// Junction counts (input level L2), when the input carries them.
    pub junctions: Option<JunctionSet>,
    /// Per-cell metadata aligned to the matrix's cell order (may be empty).
    pub metadata: CellMetadata,
}

pub fn run_stage1(
    input: &InputDescriptor,
    out_dir: &Path,
) -> Result<MmapExpressionMatrix, InputError> {
    Ok(run_stage1_full(input, out_dir, None)?.matrix)
}

/// Full stage 1: main matrix, optional layers, optional metadata.
/// `metadata_override` points at a `metadata.tsv[.gz]`; otherwise the file is
/// looked up next to a 10x directory / shared cache, or read from `obs` of an
/// AnnData input.
pub fn run_stage1_full(
    input: &InputDescriptor,
    out_dir: &Path,
    metadata_override: Option<&Path>,
) -> Result<Stage1Output, InputError> {
    let expr_path = out_dir.join("expr.bin");
    std::fs::create_dir_all(out_dir).map_err(|e| InputError::io(out_dir, e))?;

    let raw = match &input.kind {
        InputKind::TenX(tenx) => mtx::read_tenx(tenx)?,
        InputKind::H5AD(h5ad_input) => h5ad::read_h5ad(h5ad_input)?,
        InputKind::OrganelleCache(shared) => {
            info!("using shared cache mmap: {}", shared.cache_path.display());
            let matrix = MmapExpressionMatrix::open_shared_cache(&shared.cache_path)?;
            // The shared cache is already in canonical order, so layers found
            // next to it are matched by barcode against that order.
            let genes: Vec<String> =
                (0..matrix.n_genes()).map(|g| matrix.gene_symbol(g).to_string()).collect();
            let cells: Vec<String> =
                (0..matrix.n_cells()).map(|c| matrix.cell_name(c).to_string()).collect();
            let layers = match &input.layers {
                Some(location) => {
                    let raw_layers = read_layers(location, &genes, &cells)?;
                    Some(build_layers(raw_layers, genes.len(), cells.len(), None, None))
                }
                None => None,
            };
            let junctions = match &input.junctions {
                Some(location) => Some(build_junctions(read_junctions(location, &cells)?, cells.len(), None)),
                None => None,
            };
            let metadata = load_tsv_metadata(&shared.root, metadata_override, &cells)?;
            return Ok(Stage1Output {
                matrix,
                layers,
                junctions,
                metadata,
            });
        }
    };

    if raw.genes.len() != input.n_genes || raw.cells.len() != input.n_cells {
        return Err(InputError::InvalidSparseMatrix);
    }

    // Layers are read in the raw index space and reindexed with the same maps
    // as the main matrix below.
    let raw_layers = match &input.layers {
        Some(location) => Some(read_layers(location, &raw.genes, &raw.cells)?),
        None => None,
    };
    let raw_junctions = match &input.junctions {
        Some(location) => Some(read_junctions(location, &raw.cells)?),
        None => None,
    };

    let raw_cell_order = raw.cells.clone();
    let gene_index = build_index(raw.genes, true)?;
    let cell_index = build_index(raw.cells, false)?;

    let metadata = match &input.kind {
        InputKind::TenX(tenx) => {
            load_tsv_metadata(&tenx.root, metadata_override, &cell_index.sorted_names)?
        }
        InputKind::H5AD(h5) => match metadata_override {
            Some(p) => read_metadata_tsv(p, &cell_index.sorted_names)?,
            None => read_metadata_h5ad(&h5.path, &raw_cell_order, &cell_index.sorted_names)?,
        },
        InputKind::OrganelleCache(_) => unreachable!("handled above"),
    };

    let layers = raw_layers.map(|rl| {
        build_layers(
            rl,
            input.n_genes,
            input.n_cells,
            Some(&gene_index.old_to_new),
            Some(&cell_index.old_to_new),
        )
    });
    let junctions = raw_junctions.map(|rj| build_junctions(rj, input.n_cells, Some(&cell_index.old_to_new)));

    let mut triplets: Vec<(u32, u32, u32)> = raw
        .triplets
        .into_iter()
        .map(|(gene, cell, count)| {
            let new_gene = gene_index.old_to_new[gene as usize];
            let new_cell = cell_index.old_to_new[cell as usize];
            (new_gene, new_cell, count)
        })
        .collect();

    triplets.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
    let triplets = consolidate_triplets(triplets)?;

    let mut libsizes = vec![0u64; input.n_cells];
    for (_, cell, count) in &triplets {
        let cell_idx = *cell as usize;
        if cell_idx >= libsizes.len() {
            return Err(InputError::InvalidSparseMatrix);
        }
        libsizes[cell_idx] += *count as u64;
    }

    let data = CacheData {
        n_genes: input.n_genes,
        n_cells: input.n_cells,
        gene_symbols: gene_index.sorted_names,
        cell_names: cell_index.sorted_names,
        triplets,
        libsizes,
    };

    write_expr_bin(&expr_path, data)?;
    info!("expression cache written: {}", expr_path.display());

    let matrix = MmapExpressionMatrix::open(&expr_path)?;
    Ok(Stage1Output {
        matrix,
        layers,
        junctions,
        metadata,
    })
}

/// Reindexes raw junction triplets to the canonical cell order.
fn build_junctions(raw: RawJunctions, n_cells: usize, cell_map: Option<&[u32]>) -> JunctionSet {
    let n_junctions = raw.junctions.len();
    let triplets: Vec<(u32, u32, u32)> = raw
        .triplets
        .into_iter()
        .map(|(j, c, k)| (j, cell_map.map_or(c, |m| m[c as usize]), k))
        .collect();
    let counts = LayerMatrix::from_triplets(n_junctions, n_cells, triplets);
    info!(
        source = raw.source.as_str(),
        junctions = n_junctions,
        nnz = counts.nnz(),
        cells_without_junctions = raw.cells_without_junctions,
        "junction matrix indexed"
    );
    JunctionSet {
        junctions: raw.junctions,
        counts,
        source: raw.source,
        cells_without_junctions: raw.cells_without_junctions,
    }
}

/// `metadata.tsv[.gz]` next to a 10x directory (or an explicit path); empty
/// metadata when there is none.
fn load_tsv_metadata(
    root: &Path,
    override_path: Option<&Path>,
    cell_names: &[String],
) -> Result<CellMetadata, InputError> {
    if let Some(p) = override_path {
        return read_metadata_tsv(p, cell_names);
    }
    for name in ["metadata.tsv", "metadata.tsv.gz"] {
        let p = root.join(name);
        if p.is_file() {
            return read_metadata_tsv(&p, cell_names);
        }
    }
    Ok(CellMetadata::default())
}

/// Reindexes raw layer triplets (optional `old_to_new` maps) and builds the
/// CSC layer set.
fn build_layers(
    raw: RawLayers,
    n_genes: usize,
    n_cells: usize,
    gene_map: Option<&[u32]>,
    cell_map: Option<&[u32]>,
) -> SplicedUnspliced {
    let remap = |layer: RawLayer| -> LayerMatrix {
        let triplets: Vec<(u32, u32, u32)> = layer
            .triplets
            .into_iter()
            .map(|(g, c, k)| {
                let g = gene_map.map_or(g, |m| m[g as usize]);
                let c = cell_map.map_or(c, |m| m[c as usize]);
                (g, c, k)
            })
            .collect();
        LayerMatrix::from_triplets(n_genes, n_cells, triplets)
    };
    let spliced = remap(raw.spliced);
    let unspliced = remap(raw.unspliced);
    let ambiguous = raw.ambiguous.map(remap);
    info!(
        source = raw.source.as_str(),
        spliced_nnz = spliced.nnz(),
        unspliced_nnz = unspliced.nnz(),
        cells_without_layers = raw.cells_without_layers,
        "spliced/unspliced layers indexed"
    );
    SplicedUnspliced {
        spliced,
        unspliced,
        ambiguous,
        source: raw.source,
        cells_without_layers: raw.cells_without_layers,
    }
}

fn consolidate_triplets(
    mut triplets: Vec<(u32, u32, u32)>,
) -> Result<Vec<(u32, u32, u32)>, InputError> {
    if triplets.is_empty() {
        return Ok(triplets);
    }
    let mut consolidated = Vec::with_capacity(triplets.len());
    let mut current = triplets[0];
    for entry in triplets.drain(1..) {
        if entry.0 == current.0 && entry.1 == current.1 {
            let sum = current.2 as u64 + entry.2 as u64;
            if sum > u32::MAX as u64 {
                return Err(InputError::InvalidSparseMatrix);
            }
            current.2 = sum as u32;
        } else {
            consolidated.push(current);
            current = entry;
        }
    }
    consolidated.push(current);
    Ok(consolidated)
}
