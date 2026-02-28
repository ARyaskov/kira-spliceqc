use std::time::Instant;

use rayon::prelude::*;
use tracing::{debug, info, warn};

use crate::expression::ExpressionMatrix;
use crate::genesets::catalog::default_catalog_path;
use crate::genesets::{GenesetCatalog, load_catalog};
use crate::input::error::InputError;
use crate::model::geneset_activity::GenesetActivityMatrix;

pub fn run_stage2(matrix: &dyn ExpressionMatrix) -> Result<GenesetActivityMatrix, InputError> {
    let catalog_path = default_catalog_path();
    let catalog = load_catalog(&catalog_path, matrix)?;
    aggregate(matrix, &catalog)
}

pub fn aggregate(
    matrix: &dyn ExpressionMatrix,
    catalog: &GenesetCatalog,
) -> Result<GenesetActivityMatrix, InputError> {
    let n_cells = matrix.n_cells();
    let n_genesets = catalog.genesets.len();
    let mut values = vec![0.0f32; n_genesets * n_cells];
    let mut genesets = Vec::with_capacity(n_genesets);
    let mut axes = Vec::with_capacity(n_genesets);
    let libsize_scale = (0..n_cells)
        .map(|cell| 1e4_f32 / matrix.libsize(cell).max(1) as f32)
        .collect::<Vec<_>>();

    let start = Instant::now();
    for (gs_idx, geneset) in catalog.genesets.iter().enumerate() {
        genesets.push(geneset.id.clone());
        axes.push(geneset.axis.clone());
        let row_start = gs_idx * n_cells;
        let row = &mut values[row_start..row_start + n_cells];
        if geneset.gene_ids.is_empty() {
            row.fill(f32::NAN);
            warn!(
                geneset_id = geneset.id.as_str(),
                "geneset has no resolved genes"
            );
            continue;
        }

        row.par_iter_mut().enumerate().for_each(|(cell, out)| {
            let mut sum = 0.0f32;
            let scale = libsize_scale[cell];
            for &gene_id in &geneset.gene_ids {
                let count = matrix.count(gene_id as usize, cell) as f32;
                sum += (1.0 + scale * count).ln();
            }
            *out = sum / geneset.gene_ids.len() as f32;
        });

        info!(
            geneset_id = geneset.id.as_str(),
            resolved = geneset.gene_ids.len(),
            "geneset aggregated"
        );
    }

    debug!(
        elapsed_ms = start.elapsed().as_millis(),
        "geneset aggregation complete"
    );

    Ok(GenesetActivityMatrix {
        genesets,
        axes,
        values,
        n_cells,
    })
}
