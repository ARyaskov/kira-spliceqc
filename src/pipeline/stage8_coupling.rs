use std::collections::HashMap;
use std::time::Instant;

use tracing::{debug, info, warn};

use crate::input::error::InputError;
use crate::model::coupling::CouplingStressMetrics;
use crate::model::geneset_activity::GenesetActivityMatrix;
use crate::stats::robust::{mad, median};

const EPS: f32 = 1e-6;

const REQUIRED: &[&str] = &["TRANSCRIPTION_COUPLING", "U1_CORE", "U2_CORE", "SF3B_AXIS"];

pub fn run_stage8(activity: &GenesetActivityMatrix) -> Result<CouplingStressMetrics, InputError> {
    compute(activity)
}

pub fn compute(activity: &GenesetActivityMatrix) -> Result<CouplingStressMetrics, InputError> {
    info!(genesets = ?activity.genesets, "coupling axis detected");

    let n_cells = activity.n_cells;
    let mut id_to_idx = HashMap::new();
    for (idx, id) in activity.genesets.iter().enumerate() {
        id_to_idx.insert(id.as_str(), idx);
    }

    let mut z_scores: HashMap<&'static str, Vec<f32>> = HashMap::new();
    for &id in REQUIRED {
        if let Some(&idx) = id_to_idx.get(id) {
            let values = extract_geneset(activity, idx);
            let z = robust_z(&values);
            let has_finite = z.iter().any(|v| v.is_finite());
            if !has_finite {
                warn!(geneset_id = id, "geneset has no resolved genes (all NaN)");
            }
            z_scores.insert(id, z);
        } else {
            warn!(geneset_id = id, "missing geneset");
        }
    }

    let z_tx = match z_scores.remove("TRANSCRIPTION_COUPLING") {
        Some(v) if v.iter().any(|v| v.is_finite()) => v,
        _ => return Err(InputError::EmptyCouplingGeneset),
    };

    let z_u1 = z_scores
        .remove("U1_CORE")
        .unwrap_or_else(|| vec![f32::NAN; n_cells]);
    let z_u2 = z_scores
        .remove("U2_CORE")
        .unwrap_or_else(|| vec![f32::NAN; n_cells]);
    let z_sf3b = z_scores
        .remove("SF3B_AXIS")
        .unwrap_or_else(|| vec![f32::NAN; n_cells]);

    let start = Instant::now();
    let mut coupling_stress = vec![f32::NAN; n_cells];
    for cell in 0..n_cells {
        let mean_core = mean3(&z_u1, &z_u2, &z_sf3b, cell);
        let tx = z_tx[cell];
        if tx.is_finite() && mean_core.is_finite() {
            coupling_stress[cell] = tx - mean_core;
        } else {
            coupling_stress[cell] = f32::NAN;
        }
    }

    debug!(
        elapsed_ms = start.elapsed().as_millis(),
        "coupling stress computed"
    );

    Ok(CouplingStressMetrics { coupling_stress })
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

fn mean3(a: &[f32], b: &[f32], c: &[f32], cell: usize) -> f32 {
    let va = a[cell];
    let vb = b[cell];
    let vc = c[cell];
    if va.is_finite() && vb.is_finite() && vc.is_finite() {
        (va + vb + vc) / 3.0
    } else {
        f32::NAN
    }
}
