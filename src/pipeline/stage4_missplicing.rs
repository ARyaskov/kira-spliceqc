use std::time::Instant;

use ahash::AHashMap;
use rayon::prelude::*;
use tracing::{debug, info, warn};

use crate::input::error::InputError;
use crate::model::geneset_activity::GenesetActivityMatrix;
use crate::model::missplicing::MissplicingMetrics;
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

pub fn run_stage4(activity: &GenesetActivityMatrix) -> Result<MissplicingMetrics, InputError> {
    compute(activity)
}

pub fn compute(activity: &GenesetActivityMatrix) -> Result<MissplicingMetrics, InputError> {
    debug!(genesets = ?activity.genesets, "geneset ids detected");

    let n_cells = activity.n_cells;
    let id_to_idx = build_id_index(activity);

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

    if present_core < 2 {
        return Err(InputError::InsufficientCoreGenesets);
    }

    let start = Instant::now();

    let z_u1 = z_scores.get("U1_CORE");
    let z_u2 = z_scores.get("U2_CORE");
    let z_sf3b = z_scores.get("SF3B_AXIS");
    let z_srsf = z_scores.get("SRSF_SR");
    let z_hnrnp = z_scores.get("HNRNP");
    let z_u12 = z_scores.get("MINOR_U12");
    let z_nmd = z_scores.get("NMD_SURVEILLANCE");

    // Per-cell computations are independent; parallel map preserves order.
    let combined: Vec<(f32, f32, f32, f32, f32, f32)> = (0..n_cells)
        .into_par_iter()
        .map(|cell| {
            let core_mean = mean3_opt(z_u1, z_u2, z_sf3b, cell);
            let b_core = relu(-core_mean);

            let u12_val = value_at_opt(z_u12, cell);
            let major_mean = mean2_opt(z_u1, z_u2, cell);
            let b_u12 = relu(u12_val - major_mean);

            let b_nmd = relu(value_at_opt(z_nmd, cell));

            let srsf_val = value_at_opt(z_srsf, cell);
            let hnrnp_val = value_at_opt(z_hnrnp, cell);
            let b_srhn = if srsf_val.is_finite() && hnrnp_val.is_finite() {
                (srsf_val - hnrnp_val).abs()
            } else {
                f32::NAN
            };

            let components = [b_core, b_u12, b_nmd, b_srhn];
            let (burden, burden_star) = if components.iter().all(|v| v.is_finite()) {
                let mb = 0.35 * components[0]
                    + 0.25 * components[1]
                    + 0.25 * components[2]
                    + 0.15 * components[3];
                (mb, 1.0 - (-mb).exp())
            } else {
                (f32::NAN, f32::NAN)
            };

            (b_core, b_u12, b_nmd, b_srhn, burden, burden_star)
        })
        .collect();

    let mut b_core = Vec::with_capacity(n_cells);
    let mut b_u12 = Vec::with_capacity(n_cells);
    let mut b_nmd = Vec::with_capacity(n_cells);
    let mut b_srhn = Vec::with_capacity(n_cells);
    let mut burden = Vec::with_capacity(n_cells);
    let mut burden_star = Vec::with_capacity(n_cells);
    for (c, u, n, s, b, bs) in combined {
        b_core.push(c);
        b_u12.push(u);
        b_nmd.push(n);
        b_srhn.push(s);
        burden.push(b);
        burden_star.push(bs);
    }

    info!(
        elapsed_ms = start.elapsed().as_millis(),
        "missplicing metrics computed"
    );

    Ok(MissplicingMetrics {
        b_core,
        b_u12,
        b_nmd,
        b_srhn,
        burden,
        burden_star,
    })
}

pub(crate) fn build_id_index(activity: &GenesetActivityMatrix) -> AHashMap<&str, usize> {
    let mut map = AHashMap::with_capacity(activity.genesets.len());
    for (idx, id) in activity.genesets.iter().enumerate() {
        map.insert(id.as_str(), idx);
    }
    map
}

#[inline]
fn value_at_opt(values: Option<&Vec<f32>>, cell: usize) -> f32 {
    match values {
        Some(v) => v[cell],
        None => f32::NAN,
    }
}

#[inline]
fn mean2_opt(a: Option<&Vec<f32>>, b: Option<&Vec<f32>>, cell: usize) -> f32 {
    let va = value_at_opt(a, cell);
    let vb = value_at_opt(b, cell);
    if va.is_finite() && vb.is_finite() {
        (va + vb) * 0.5
    } else {
        f32::NAN
    }
}

#[inline]
fn mean3_opt(a: Option<&Vec<f32>>, b: Option<&Vec<f32>>, c: Option<&Vec<f32>>, cell: usize) -> f32 {
    let va = value_at_opt(a, cell);
    let vb = value_at_opt(b, cell);
    let vc = value_at_opt(c, cell);
    if va.is_finite() && vb.is_finite() && vc.is_finite() {
        (va + vb + vc) / 3.0
    } else {
        f32::NAN
    }
}

#[inline]
fn relu(value: f32) -> f32 {
    if !value.is_finite() {
        f32::NAN
    } else if value > 0.0 {
        value
    } else {
        0.0
    }
}
