use std::time::Instant;

use ahash::AHashMap;
use rayon::prelude::*;
use tracing::{debug, info, warn};

use crate::input::error::InputError;
use crate::model::geneset_activity::GenesetActivityMatrix;
use crate::model::imbalance::SpliceosomeImbalanceMetrics;
use crate::stats::robust::{extract_geneset_slice, robust_z};

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
    debug!(genesets = ?activity.genesets, "geneset ids detected");

    let n_cells = activity.n_cells;
    let id_to_idx = super::stage4_missplicing::build_id_index(activity);

    let mut z_scores: AHashMap<&'static str, Vec<f32>> = AHashMap::with_capacity(REQUIRED_GENESETS.len());
    let mut present_core = 0usize;

    for &id in REQUIRED_GENESETS {
        if let Some(&idx) = id_to_idx.get(id) {
            let slice = extract_geneset_slice(&activity.values, idx, n_cells);
            let (z, _) = robust_z(slice);
            let has_finite = z.iter().any(|v| v.is_finite());
            if !has_finite {
                warn!(geneset_id = id, "geneset has no resolved genes (all NaN)");
            }
            if has_finite && matches!(id, "U1_CORE" | "U2_CORE" | "SF3B_AXIS") {
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

    let z_u1 = z_scores.remove("U1_CORE").unwrap_or_else(|| vec![f32::NAN; n_cells]);
    let z_u2 = z_scores.remove("U2_CORE").unwrap_or_else(|| vec![f32::NAN; n_cells]);
    let z_sf3b = z_scores.remove("SF3B_AXIS").unwrap_or_else(|| vec![f32::NAN; n_cells]);
    let z_srsf = z_scores.remove("SRSF_SR").unwrap_or_else(|| vec![f32::NAN; n_cells]);
    let z_hnrnp = z_scores.remove("HNRNP").unwrap_or_else(|| vec![f32::NAN; n_cells]);
    let z_u12 = z_scores.remove("MINOR_U12").unwrap_or_else(|| vec![f32::NAN; n_cells]);
    let z_nmd = z_scores.remove("NMD_SURVEILLANCE").unwrap_or_else(|| vec![f32::NAN; n_cells]);

    let derived: Vec<(f32, f32, f32, f32, f32)> = (0..n_cells)
        .into_par_iter()
        .map(|cell| {
            let u1 = z_u1[cell];
            let u2 = z_u2[cell];
            let sf3b = z_sf3b[cell];
            let srsf = z_srsf[cell];
            let hnrnp = z_hnrnp[cell];
            let u12 = z_u12[cell];
            let nmd = z_nmd[cell];

            let sr_hnrnp = if srsf.is_finite() && hnrnp.is_finite() {
                srsf - hnrnp
            } else {
                f32::NAN
            };
            let u2_u1 = if u2.is_finite() && u1.is_finite() {
                u2 - u1
            } else {
                f32::NAN
            };
            let u12_major = if u12.is_finite() && u1.is_finite() && u2.is_finite() {
                u12 - (u1 + u2) * 0.5
            } else {
                f32::NAN
            };
            let nmd_axis = if nmd.is_finite() { nmd } else { f32::NAN };

            let mut vals = [u1, u2, sf3b, srsf, hnrnp, u12];
            let imbalance = if vals.iter().all(|v| v.is_finite()) {
                for v in &mut vals {
                    *v = v.clamp(-6.0, 6.0);
                }
                let mean_sq: f32 = vals.iter().map(|v| v * v).sum::<f32>() / vals.len() as f32;
                mean_sq.sqrt()
            } else {
                f32::NAN
            };

            (sr_hnrnp, u2_u1, u12_major, nmd_axis, imbalance)
        })
        .collect();

    let (axis_sr_hnrnp, axis_u2_u1, axis_u12_major, axis_nmd, imbalance) = unzip5(derived);

    info!(
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

#[inline]
fn unzip5<A, B, C, D, E>(
    src: Vec<(A, B, C, D, E)>,
) -> (Vec<A>, Vec<B>, Vec<C>, Vec<D>, Vec<E>) {
    let n = src.len();
    let mut a = Vec::with_capacity(n);
    let mut b = Vec::with_capacity(n);
    let mut c = Vec::with_capacity(n);
    let mut d = Vec::with_capacity(n);
    let mut e = Vec::with_capacity(n);
    for (av, bv, cv, dv, ev) in src {
        a.push(av);
        b.push(bv);
        c.push(cv);
        d.push(dv);
        e.push(ev);
    }
    (a, b, c, d, e)
}
