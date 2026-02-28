use std::path::Path;

use tracing::info;

use crate::expression::cache_writer::{CacheData, write_expr_bin};
use crate::expression::index::build_index;
use crate::expression::mmap::MmapExpressionMatrix;
use crate::input::error::InputError;
use crate::input::{InputDescriptor, InputKind};
use crate::io::{h5ad, mtx};

pub fn run_stage1(
    input: &InputDescriptor,
    out_dir: &Path,
) -> Result<MmapExpressionMatrix, InputError> {
    let expr_path = out_dir.join("expr.bin");
    std::fs::create_dir_all(out_dir).map_err(|e| InputError::io(out_dir, e))?;

    let raw = match &input.kind {
        InputKind::TenX(tenx) => mtx::read_tenx(tenx)?,
        InputKind::H5AD(h5ad_input) => h5ad::read_h5ad(h5ad_input)?,
        InputKind::OrganelleCache(shared) => {
            info!("using shared cache mmap: {}", shared.cache_path.display());
            return MmapExpressionMatrix::open_shared_cache(&shared.cache_path);
        }
    };

    if raw.genes.len() != input.n_genes || raw.cells.len() != input.n_cells {
        return Err(InputError::InvalidSparseMatrix);
    }

    let gene_index = build_index(raw.genes, true)?;
    let cell_index = build_index(raw.cells, false)?;

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

    MmapExpressionMatrix::open(&expr_path)
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
