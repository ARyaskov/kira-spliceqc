use std::collections::HashMap;
use std::time::Instant;

use tracing::{debug, info, warn};

use crate::input::error::InputError;
use crate::model::geneset_activity::GenesetActivityMatrix;
use crate::model::missplicing::MissplicingMetrics;
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

pub fn run_stage4(activity: &GenesetActivityMatrix) -> Result<MissplicingMetrics, InputError> {
    compute(activity)
}

pub fn compute(activity: &GenesetActivityMatrix) -> Result<MissplicingMetrics, InputError> {
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
            if matches!(id, "U1_CORE" | "U2_CORE" | "SF3B_AXIS") {
                // no increment
            }
        }
    }

    if present_core < 2 {
        return Err(InputError::InsufficientCoreGenesets);
    }

    let start = Instant::now();

    let mut b_core = vec![f32::NAN; n_cells];
    let mut b_u12 = vec![f32::NAN; n_cells];
    let mut b_nmd = vec![f32::NAN; n_cells];
    let mut b_srhn = vec![f32::NAN; n_cells];
    let mut burden = vec![f32::NAN; n_cells];
    let mut burden_star = vec![f32::NAN; n_cells];

    let z_u1 = z_scores.get("U1_CORE");
    let z_u2 = z_scores.get("U2_CORE");
    let z_sf3b = z_scores.get("SF3B_AXIS");
    let z_srsf = z_scores.get("SRSF_SR");
    let z_hnrnp = z_scores.get("HNRNP");
    let z_u12 = z_scores.get("MINOR_U12");
    let z_nmd = z_scores.get("NMD_SURVEILLANCE");

    for cell in 0..n_cells {
        let core_mean = mean3(z_u1, z_u2, z_sf3b, cell);
        b_core[cell] = relu(-core_mean);

        let u12_val = value_at(z_u12, cell);
        let major_mean = mean2(z_u1, z_u2, cell);
        b_u12[cell] = relu(u12_val - major_mean);

        b_nmd[cell] = relu(value_at(z_nmd, cell));

        let srsf_val = value_at(z_srsf, cell);
        let hnrnp_val = value_at(z_hnrnp, cell);
        b_srhn[cell] = if srsf_val.is_finite() && hnrnp_val.is_finite() {
            (srsf_val - hnrnp_val).abs()
        } else {
            f32::NAN
        };

        let components = [b_core[cell], b_u12[cell], b_nmd[cell], b_srhn[cell]];
        if components.iter().all(|v| v.is_finite()) {
            let mb = 0.35 * components[0]
                + 0.25 * components[1]
                + 0.25 * components[2]
                + 0.15 * components[3];
            burden[cell] = mb;
            burden_star[cell] = 1.0 - (-mb).exp();
        } else {
            burden[cell] = f32::NAN;
            burden_star[cell] = f32::NAN;
        }
    }

    debug!(
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

fn value_at(values: Option<&Vec<f32>>, cell: usize) -> f32 {
    match values {
        Some(v) => v[cell],
        None => f32::NAN,
    }
}

fn mean2(a: Option<&Vec<f32>>, b: Option<&Vec<f32>>, cell: usize) -> f32 {
    let va = value_at(a, cell);
    let vb = value_at(b, cell);
    if va.is_finite() && vb.is_finite() {
        (va + vb) * 0.5
    } else {
        f32::NAN
    }
}

fn mean3(a: Option<&Vec<f32>>, b: Option<&Vec<f32>>, c: Option<&Vec<f32>>, cell: usize) -> f32 {
    let va = value_at(a, cell);
    let vb = value_at(b, cell);
    let vc = value_at(c, cell);
    if va.is_finite() && vb.is_finite() && vc.is_finite() {
        (va + vb + vc) / 3.0
    } else {
        f32::NAN
    }
}

fn relu(value: f32) -> f32 {
    if !value.is_finite() {
        f32::NAN
    } else if value > 0.0 {
        value
    } else {
        0.0
    }
}
