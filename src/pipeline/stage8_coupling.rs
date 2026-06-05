use std::time::Instant;

use ahash::AHashMap;
use rayon::prelude::*;
use tracing::{debug, info, warn};

use crate::input::error::InputError;
use crate::model::coupling::CouplingStressMetrics;
use crate::model::geneset_activity::GenesetActivityMatrix;
use crate::stats::robust::{extract_geneset_slice, robust_z};

const REQUIRED: &[&str] = &["TRANSCRIPTION_COUPLING", "U1_CORE", "U2_CORE", "SF3B_AXIS"];

pub fn run_stage8(activity: &GenesetActivityMatrix) -> Result<CouplingStressMetrics, InputError> {
    compute(activity)
}

pub fn compute(activity: &GenesetActivityMatrix) -> Result<CouplingStressMetrics, InputError> {
    debug!(genesets = ?activity.genesets, "coupling axis detected");

    let n_cells = activity.n_cells;
    let id_to_idx = super::stage4_missplicing::build_id_index(activity);

    let mut z_scores: AHashMap<&'static str, Vec<f32>> = AHashMap::with_capacity(REQUIRED.len());
    for &id in REQUIRED {
        if let Some(&idx) = id_to_idx.get(id) {
            let slice = extract_geneset_slice(&activity.values, idx, n_cells);
            let (z, _) = robust_z(slice);
            if !z.iter().any(|v| v.is_finite()) {
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

    let z_u1 = z_scores.remove("U1_CORE").unwrap_or_else(|| vec![f32::NAN; n_cells]);
    let z_u2 = z_scores.remove("U2_CORE").unwrap_or_else(|| vec![f32::NAN; n_cells]);
    let z_sf3b = z_scores.remove("SF3B_AXIS").unwrap_or_else(|| vec![f32::NAN; n_cells]);

    let start = Instant::now();
    let coupling_stress: Vec<f32> = (0..n_cells)
        .into_par_iter()
        .map(|cell| {
            let core = mean3(&z_u1, &z_u2, &z_sf3b, cell);
            let tx = z_tx[cell];
            if tx.is_finite() && core.is_finite() {
                tx - core
            } else {
                f32::NAN
            }
        })
        .collect();

    info!(
        elapsed_ms = start.elapsed().as_millis(),
        "coupling stress computed"
    );

    Ok(CouplingStressMetrics { coupling_stress })
}

#[inline]
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
