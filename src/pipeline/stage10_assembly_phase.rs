use std::collections::HashMap;
use std::time::Instant;

use tracing::{debug, info, warn};

use crate::input::error::InputError;
use crate::model::assembly_phase::AssemblyPhaseImbalanceMetrics;
use crate::model::geneset_activity::GenesetActivityMatrix;
use crate::stats::robust::{mad, median};

const EPS: f32 = 1e-6;

const REQUIRED: &[&str] = &[
    "SPLICE_EA_PHASE",
    "SPLICE_B_PHASE",
    "SPLICE_CATALYTIC_PHASE",
];

pub fn run_stage10(
    activity: &GenesetActivityMatrix,
) -> Result<AssemblyPhaseImbalanceMetrics, InputError> {
    compute(activity)
}

pub fn compute(
    activity: &GenesetActivityMatrix,
) -> Result<AssemblyPhaseImbalanceMetrics, InputError> {
    info!(genesets = ?activity.genesets, "assembly phase genesets resolved");

    let n_cells = activity.n_cells;
    let mut id_to_idx = HashMap::new();
    for (idx, id) in activity.genesets.iter().enumerate() {
        id_to_idx.insert(id.as_str(), idx);
    }

    let mut z_scores: HashMap<&'static str, Vec<f32>> = HashMap::new();
    let mut present = 0usize;

    for &id in REQUIRED {
        if let Some(&idx) = id_to_idx.get(id) {
            let values = extract_geneset(activity, idx);
            let z = robust_z(&values);
            let has_finite = z.iter().any(|v| v.is_finite());
            if !has_finite {
                warn!(geneset_id = id, "geneset has no resolved genes (all NaN)");
            }
            if has_finite {
                present += 1;
            }
            z_scores.insert(id, z);
        } else {
            warn!(geneset_id = id, "missing geneset");
        }
    }

    if present < 2 {
        return Err(InputError::InsufficientAssemblyPhaseGenesets);
    }

    let start = Instant::now();

    let z_ea = z_scores
        .remove("SPLICE_EA_PHASE")
        .unwrap_or_else(|| vec![f32::NAN; n_cells]);
    let z_b = z_scores
        .remove("SPLICE_B_PHASE")
        .unwrap_or_else(|| vec![f32::NAN; n_cells]);
    let z_cat = z_scores
        .remove("SPLICE_CATALYTIC_PHASE")
        .unwrap_or_else(|| vec![f32::NAN; n_cells]);

    let mut ea_imbalance = vec![f32::NAN; n_cells];
    let mut b_imbalance = vec![f32::NAN; n_cells];
    let mut cat_imbalance = vec![f32::NAN; n_cells];

    for cell in 0..n_cells {
        let ea = z_ea[cell];
        let b = z_b[cell];
        let cat = z_cat[cell];
        if ea.is_finite() && b.is_finite() && cat.is_finite() {
            ea_imbalance[cell] = ea - (b + cat) * 0.5;
            b_imbalance[cell] = b - (ea + cat) * 0.5;
            cat_imbalance[cell] = cat - (ea + b) * 0.5;
        } else {
            ea_imbalance[cell] = f32::NAN;
            b_imbalance[cell] = f32::NAN;
            cat_imbalance[cell] = f32::NAN;
        }
    }

    debug!(
        elapsed_ms = start.elapsed().as_millis(),
        "assembly phase imbalance computed"
    );

    Ok(AssemblyPhaseImbalanceMetrics {
        z_ea,
        z_b,
        z_cat,
        ea_imbalance,
        b_imbalance,
        cat_imbalance,
    })
}

fn extract_geneset(activity: &GenesetActivityMatrix, idx: usize) -> Vec<f32> {
    let n_cells = activity.n_cells;
    let mut values = Vec::with_capacity(n_cells);
    let start = idx * n_cells;
    let end = start + n_cells;
    values.extend_from_slice(&activity.values[start..end]);
    values
}

fn robust_z(values: &[f32]) -> Vec<f32> {
    let med = median(values);
    let mad_val = mad(values, med) * 1.4826 + EPS;
    values
        .iter()
        .map(|&v| {
            if v.is_finite() && med.is_finite() && mad_val.is_finite() {
                (v - med) / mad_val
            } else {
                f32::NAN
            }
        })
        .collect()
}
