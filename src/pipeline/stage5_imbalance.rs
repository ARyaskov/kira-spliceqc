use std::collections::HashMap;
use std::time::Instant;

use tracing::{debug, info, warn};

use crate::input::error::InputError;
use crate::model::geneset_activity::GenesetActivityMatrix;
use crate::model::imbalance::SpliceosomeImbalanceMetrics;
use crate::stats::robust::{mad, median};

const EPS: f32 = 1e-6;

const REQUIRED_GENESETS: &[&str] = &[
    "U1_CORE",
    "U2_CORE",
    "SF3B_AXIS",
    "SRSF_SR",
    "HNRNP",
    "MINOR_U12",
    "NMD_SURVEILLANCE",
];

pub fn run_stage5(
    activity: &GenesetActivityMatrix,
) -> Result<SpliceosomeImbalanceMetrics, InputError> {
    compute(activity)
}

pub fn compute(
    activity: &GenesetActivityMatrix,
) -> Result<SpliceosomeImbalanceMetrics, InputError> {
    info!(genesets = ?activity.genesets, "geneset ids detected");

    let n_cells = activity.n_cells;
    let mut id_to_idx = HashMap::new();
    for (idx, id) in activity.genesets.iter().enumerate() {
        id_to_idx.insert(id.as_str(), idx);
    }

    let mut z_scores: HashMap<&'static str, Vec<f32>> = HashMap::new();
    let mut present_core = 0usize;

    for &id in REQUIRED_GENESETS {
        if let Some(&idx) = id_to_idx.get(id) {
            let values = extract_geneset(activity, idx);
            let z = robust_z(&values);
            let has_finite = z.iter().any(|v| v.is_finite());
            if !has_finite {
                warn!(geneset_id = id, "geneset has no resolved genes (all NaN)");
            }
            if matches!(id, "U1_CORE" | "U2_CORE" | "SF3B_AXIS") && has_finite {
                present_core += 1;
            }
            z_scores.insert(id, z);
        } else {
            warn!(geneset_id = id, "missing geneset");
        }
    }

    if present_core < 3 {
        return Err(InputError::InsufficientCoreGenesetsImbalance);
    }

    let start = Instant::now();

    let z_u1 = z_scores
        .remove("U1_CORE")
        .unwrap_or_else(|| vec![f32::NAN; n_cells]);
    let z_u2 = z_scores
        .remove("U2_CORE")
        .unwrap_or_else(|| vec![f32::NAN; n_cells]);
    let z_sf3b = z_scores
        .remove("SF3B_AXIS")
        .unwrap_or_else(|| vec![f32::NAN; n_cells]);
    let z_srsf = z_scores
        .remove("SRSF_SR")
        .unwrap_or_else(|| vec![f32::NAN; n_cells]);
    let z_hnrnp = z_scores
        .remove("HNRNP")
        .unwrap_or_else(|| vec![f32::NAN; n_cells]);
    let z_u12 = z_scores
        .remove("MINOR_U12")
        .unwrap_or_else(|| vec![f32::NAN; n_cells]);
    let z_nmd = z_scores
        .remove("NMD_SURVEILLANCE")
        .unwrap_or_else(|| vec![f32::NAN; n_cells]);

    let mut axis_sr_hnrnp = vec![f32::NAN; n_cells];
    let mut axis_u2_u1 = vec![f32::NAN; n_cells];
    let mut axis_u12_major = vec![f32::NAN; n_cells];
    let mut axis_nmd = vec![f32::NAN; n_cells];
    let mut imbalance = vec![f32::NAN; n_cells];

    for cell in 0..n_cells {
        let u1 = z_u1[cell];
        let u2 = z_u2[cell];
        let sf3b = z_sf3b[cell];
        let srsf = z_srsf[cell];
        let hnrnp = z_hnrnp[cell];
        let u12 = z_u12[cell];
        let nmd = z_nmd[cell];

        axis_sr_hnrnp[cell] = if srsf.is_finite() && hnrnp.is_finite() {
            srsf - hnrnp
        } else {
            f32::NAN
        };

        axis_u2_u1[cell] = if u2.is_finite() && u1.is_finite() {
            u2 - u1
        } else {
            f32::NAN
        };

        axis_u12_major[cell] = if u12.is_finite() && u1.is_finite() && u2.is_finite() {
            u12 - (u1 + u2) * 0.5
        } else {
            f32::NAN
        };

        axis_nmd[cell] = if nmd.is_finite() { nmd } else { f32::NAN };

        let mut vals = [u1, u2, sf3b, srsf, hnrnp, u12];
        if vals.iter().all(|v| v.is_finite()) {
            for v in &mut vals {
                *v = clamp(*v, -6.0, 6.0);
            }
            let mean_sq = vals.iter().map(|v| v * v).sum::<f32>() / vals.len() as f32;
            imbalance[cell] = mean_sq.sqrt();
        } else {
            imbalance[cell] = f32::NAN;
        }
    }

    debug!(
        elapsed_ms = start.elapsed().as_millis(),
        "imbalance metrics computed"
    );

    Ok(SpliceosomeImbalanceMetrics {
        z_u1,
        z_u2,
        z_sf3b,
        z_srsf,
        z_hnrnp,
        z_u12,
        z_nmd,
        axis_sr_hnrnp,
        axis_u2_u1,
        axis_u12_major,
        axis_nmd,
        imbalance,
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

fn clamp(value: f32, lo: f32, hi: f32) -> f32 {
    if value < lo {
        lo
    } else if value > hi {
        hi
    } else {
        value
    }
}
