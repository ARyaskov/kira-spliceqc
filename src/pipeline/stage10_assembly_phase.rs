use std::time::Instant;

use ahash::AHashMap;
use rayon::prelude::*;
use tracing::{debug, info, warn};

use crate::input::error::InputError;
use crate::model::assembly_phase::AssemblyPhaseImbalanceMetrics;
use crate::model::geneset_activity::GenesetActivityMatrix;
use crate::stats::robust::{extract_geneset_slice, robust_z};

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
    debug!(genesets = ?activity.genesets, "assembly phase genesets resolved");

    let n_cells = activity.n_cells;
    let id_to_idx = super::stage4_missplicing::build_id_index(activity);

    let mut z_scores: AHashMap<&'static str, Vec<f32>> = AHashMap::with_capacity(REQUIRED.len());
    let mut present = 0usize;
    for &id in REQUIRED {
        if let Some(&idx) = id_to_idx.get(id) {
            let slice = extract_geneset_slice(&activity.values, idx, n_cells);
            let (z, _) = robust_z(slice);
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

    let z_ea = z_scores.remove("SPLICE_EA_PHASE").unwrap_or_else(|| vec![f32::NAN; n_cells]);
    let z_b = z_scores.remove("SPLICE_B_PHASE").unwrap_or_else(|| vec![f32::NAN; n_cells]);
    let z_cat = z_scores.remove("SPLICE_CATALYTIC_PHASE").unwrap_or_else(|| vec![f32::NAN; n_cells]);

    let derived: Vec<(f32, f32, f32)> = (0..n_cells)
        .into_par_iter()
        .map(|cell| {
            let ea = z_ea[cell];
            let b = z_b[cell];
            let cat = z_cat[cell];
            if ea.is_finite() && b.is_finite() && cat.is_finite() {
                let ea_imb = ea - (b + cat) * 0.5;
                let b_imb = b - (ea + cat) * 0.5;
                let cat_imb = cat - (ea + b) * 0.5;
                (ea_imb, b_imb, cat_imb)
            } else {
                (f32::NAN, f32::NAN, f32::NAN)
            }
        })
        .collect();

    let mut ea_imbalance = Vec::with_capacity(n_cells);
    let mut b_imbalance = Vec::with_capacity(n_cells);
    let mut cat_imbalance = Vec::with_capacity(n_cells);
    for (a, b, c) in derived {
        ea_imbalance.push(a);
        b_imbalance.push(b);
        cat_imbalance.push(c);
    }

    info!(
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
