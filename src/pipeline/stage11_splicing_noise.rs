use std::collections::HashMap;
use std::time::Instant;

use tracing::{debug, info, warn};

use crate::input::error::InputError;
use crate::model::geneset_activity::GenesetActivityMatrix;
use crate::model::splicing_noise::SplicingNoiseMetrics;
use crate::stats::robust::{mad, median};

const EPS: f32 = 1e-6;

const CORE: &[&str] = &[
    "U1_CORE",
    "U2_CORE",
    "SF3B_AXIS",
    "SRSF_SR",
    "HNRNP",
    "MINOR_U12",
];

pub fn run_stage11(activity: &GenesetActivityMatrix) -> Result<SplicingNoiseMetrics, InputError> {
    compute(activity)
}

pub fn compute(activity: &GenesetActivityMatrix) -> Result<SplicingNoiseMetrics, InputError> {
    info!(genesets = ?activity.genesets, "splicing noise genesets resolved");

    let start = Instant::now();

    let mut id_to_idx = HashMap::new();
    for (idx, id) in activity.genesets.iter().enumerate() {
        id_to_idx.insert(id.as_str(), idx);
    }

    let mut present = 0usize;
    let mut used = 0usize;
    let mut sum = 0.0_f32;
    let mut per_geneset_noise: Vec<(String, f32)> = Vec::new();

    for &id in CORE {
        if let Some(&idx) = id_to_idx.get(id) {
            present += 1;
            let values = extract_geneset(activity, idx);
            let noise = compute_noise(id, &values);
            if noise.is_finite() {
                sum += noise;
                used += 1;
            }
            per_geneset_noise.push((id.to_string(), noise));
        } else {
            warn!(geneset_id = id, "missing geneset");
        }
    }

    if present < 2 {
        return Err(InputError::InsufficientSplicingNoiseGenesets);
    }

    let noise_index = if used > 0 {
        sum / used as f32
    } else {
        f32::NAN
    };

    per_geneset_noise.sort_by(|a, b| a.0.cmp(&b.0));

    info!(genesets_used = used, "splicing noise index computed");
    debug!(
        elapsed_ms = start.elapsed().as_millis(),
        "splicing noise computed"
    );

    Ok(SplicingNoiseMetrics {
        noise_index,
        per_geneset_noise,
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

fn compute_noise(geneset_id: &str, values: &[f32]) -> f32 {
    let has_finite = values.iter().any(|v| v.is_finite());
    if !has_finite {
        warn!(
            geneset_id = geneset_id,
            "geneset has no resolved genes (all NaN)"
        );
        return f32::NAN;
    }

    let med = median(values);
    if !med.is_finite() {
        warn!(geneset_id = geneset_id, "geneset median undefined");
        return f32::NAN;
    }
    if med.abs() <= EPS {
        warn!(geneset_id = geneset_id, "geneset median near zero");
        return f32::NAN;
    }

    let mad_val = mad(values, med);
    if !mad_val.is_finite() {
        warn!(geneset_id = geneset_id, "geneset MAD undefined");
        return f32::NAN;
    }

    mad_val / (med.abs() + EPS)
}
