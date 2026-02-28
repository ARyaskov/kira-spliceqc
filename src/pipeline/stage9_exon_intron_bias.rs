use std::collections::HashMap;
use std::time::Instant;

use tracing::{debug, info, warn};

use crate::input::error::InputError;
use crate::model::exon_intron_bias::ExonIntronDefinitionMetrics;
use crate::model::geneset_activity::GenesetActivityMatrix;
use crate::stats::robust::{mad, median};

const EPS: f32 = 1e-6;

const REQUIRED: &[&str] = &["SRSF_SR", "HNRNP", "U2AF_AXIS"];

pub fn run_stage9(
    activity: &GenesetActivityMatrix,
) -> Result<ExonIntronDefinitionMetrics, InputError> {
    compute(activity)
}

pub fn compute(
    activity: &GenesetActivityMatrix,
) -> Result<ExonIntronDefinitionMetrics, InputError> {
    info!(genesets = ?activity.genesets, "exon/intron genesets resolved");

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
        return Err(InputError::InsufficientGenesetsExonIntron);
    }

    let start = Instant::now();

    let z_srsf = z_scores
        .remove("SRSF_SR")
        .unwrap_or_else(|| vec![f32::NAN; n_cells]);
    let z_hnrnp = z_scores
        .remove("HNRNP")
        .unwrap_or_else(|| vec![f32::NAN; n_cells]);
    let z_u2af = z_scores
        .remove("U2AF_AXIS")
        .unwrap_or_else(|| vec![f32::NAN; n_cells]);

    let mut exon_definition_bias = vec![f32::NAN; n_cells];
    for cell in 0..n_cells {
        let srsf = z_srsf[cell];
        let u2af = z_u2af[cell];
        let hnrnp = z_hnrnp[cell];
        if srsf.is_finite() && u2af.is_finite() && hnrnp.is_finite() {
            exon_definition_bias[cell] = (srsf + u2af) - hnrnp;
        } else {
            exon_definition_bias[cell] = f32::NAN;
        }
    }

    debug!(
        elapsed_ms = start.elapsed().as_millis(),
        "exon/intron bias computed"
    );

    Ok(ExonIntronDefinitionMetrics {
        exon_definition_bias,
        z_srsf,
        z_u2af,
        z_hnrnp,
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
