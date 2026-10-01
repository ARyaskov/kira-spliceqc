use std::time::Instant;

use rayon::prelude::*;
use tracing::{debug, info, warn};

use crate::expression::ExpressionMatrix;
use crate::genesets::catalog::default_catalog_path;
use crate::genesets::controls::ControlPool;
use crate::genesets::{GenesetCatalog, load_catalog};
use crate::input::error::InputError;
use crate::model::geneset_activity::GenesetActivityMatrix;
use crate::reference::external::{ExpressionNormsBundle, ReferenceFile};
use crate::reference::{
    Strata, apply_depth_bin_norms, depth_bin_norms, robust_z_by_stratum,
    robust_z_by_stratum_and_depth,
};

/// Depth-adaptive, stratified standardization of a raw activity matrix:
/// every panel score becomes a robust z-score against the cells of the same
/// stratum and library-size bin (`reference::robust_z_by_stratum_and_depth`).
/// Stages 4, 5, 8, 9, 10 and 13 consume the standardized matrix; stage 11
/// (noise) keeps the raw one. Panels whose MAD collapses to zero in a bin
/// are logged and yield NaN there.
pub fn standardize_activity(
    activity: &GenesetActivityMatrix,
    strata: &Strata,
    libsize: Option<&[u64]>,
) -> GenesetActivityMatrix {
    standardize_activity_with(activity, strata, libsize, None)
}

/// `standardize_activity` taking the norms from an external reference for
/// every geneset the file knows (others fall back to the dataset's own).
pub fn standardize_activity_with(
    activity: &GenesetActivityMatrix,
    strata: &Strata,
    libsize: Option<&[u64]>,
    external: Option<&ReferenceFile>,
) -> GenesetActivityMatrix {
    let n_cells = activity.n_cells;
    let mut values = vec![f32::NAN; activity.values.len()];
    for (idx, id) in activity.genesets.iter().enumerate() {
        let slice = &activity.values[idx * n_cells..(idx + 1) * n_cells];
        let external_norms = external.and_then(|f| {
            let norms = f.expression_norms(id);
            norms.iter().any(|n| n.is_some()).then_some(norms)
        });
        let z = match (external_norms, libsize) {
            (Some(norms), Some(l)) => apply_depth_bin_norms(slice, &strata.labels, l, &norms),
            (Some(_), None) | (None, None) => robust_z_by_stratum(slice, strata).0,
            (None, Some(l)) => robust_z_by_stratum_and_depth(slice, strata, l).0,
        };
        let finite_in = slice.iter().filter(|v| v.is_finite()).count();
        let finite_out = z.iter().filter(|v| v.is_finite()).count();
        if finite_in > 0 && finite_out < finite_in {
            warn!(
                geneset_id = id.as_str(),
                undefined_after_standardization = finite_in - finite_out,
                "panel MAD is zero in some stratum/depth bin; z-scores undefined there"
            );
        }
        values[idx * n_cells..(idx + 1) * n_cells].copy_from_slice(&z);
    }
    GenesetActivityMatrix {
        genesets: activity.genesets.clone(),
        axes: activity.axes.clone(),
        values,
        n_cells,
    }
}

/// Per-stratum depth-binned norms of a raw activity matrix keyed by geneset
/// id (what `reference build` stores).
pub fn activity_norms(
    activity: &GenesetActivityMatrix,
    strata: &Strata,
    libsize: &[u64],
) -> ExpressionNormsBundle {
    let n_cells = activity.n_cells;
    activity
        .genesets
        .iter()
        .enumerate()
        .map(|(idx, id)| {
            let slice = &activity.values[idx * n_cells..(idx + 1) * n_cells];
            (id.clone(), depth_bin_norms(slice, strata, libsize))
        })
        .filter(|(_, norms)| norms.iter().any(|n| !n.norms.is_empty()))
        .collect()
}

pub fn run_stage2(matrix: &dyn ExpressionMatrix) -> Result<GenesetActivityMatrix, InputError> {
    let catalog_path = default_catalog_path();
    let catalog = load_catalog(&catalog_path, matrix)?;
    aggregate(matrix, &catalog)
}

/// Raw panel means (no background subtraction).
pub fn aggregate(
    matrix: &dyn ExpressionMatrix,
    catalog: &GenesetCatalog,
) -> Result<GenesetActivityMatrix, InputError> {
    aggregate_with_controls(matrix, catalog, None)
}

/// Panel means with the control-gene background subtracted when a
/// `ControlPool` is given: `score = mean_panel(log1p cp10k) - mean_controls(log1p cp10k)`
/// (Tirosh et al. 2016; Seurat AddModuleScore). A panel whose control set
/// is empty falls back to the raw mean with a warning.
pub fn aggregate_with_controls(
    matrix: &dyn ExpressionMatrix,
    catalog: &GenesetCatalog,
    controls: Option<&ControlPool>,
) -> Result<GenesetActivityMatrix, InputError> {
    let n_cells = matrix.n_cells();
    let n_genesets = catalog.genesets.len();
    let mut values = vec![0.0f32; n_genesets * n_cells];

    let genesets: Vec<String> = catalog.genesets.iter().map(|g| g.id.clone()).collect();
    let axes: Vec<String> = catalog.genesets.iter().map(|g| g.axis.clone()).collect();

    let libsize_scale: Vec<f32> = (0..n_cells)
        .map(|cell| 1e4_f32 / matrix.libsize(cell).max(1) as f32)
        .collect();

    let start = Instant::now();
    for (gs_idx, geneset) in catalog.genesets.iter().enumerate() {
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

        let panel = geneset.gene_ids.as_slice();
        let n_panel = panel.len() as f32;
        let control_set: Vec<u32> = controls.map(|p| p.controls_for(panel)).unwrap_or_default();
        if controls.is_some() && control_set.is_empty() {
            warn!(
                geneset_id = geneset.id.as_str(),
                "no eligible control genes; panel score is the raw mean"
            );
        }
        let n_ctrl = control_set.len() as f32;
        row.par_iter_mut().enumerate().for_each(|(cell, out)| {
            let scale = libsize_scale[cell];
            let sum = matrix.panel_ln1p_scaled_sum(panel, cell, scale);
            let background = if control_set.is_empty() {
                0.0
            } else {
                matrix.panel_ln1p_scaled_sum(&control_set, cell, scale) / n_ctrl
            };
            *out = sum / n_panel - background;
        });

        debug!(
            geneset_id = geneset.id.as_str(),
            resolved = geneset.gene_ids.len(),
            controls = control_set.len(),
            "geneset aggregated"
        );
    }

    info!(
        elapsed_ms = start.elapsed().as_millis(),
        genesets = n_genesets,
        "geneset aggregation complete"
    );

    Ok(GenesetActivityMatrix {
        genesets,
        axes,
        values,
        n_cells,
    })
}
