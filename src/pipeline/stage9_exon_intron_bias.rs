use std::time::Instant;

use ahash::AHashMap;
use rayon::prelude::*;
use tracing::{debug, info, warn};

use crate::input::error::InputError;
use crate::model::exon_intron_bias::ExonIntronDefinitionMetrics;
use crate::model::geneset_activity::GenesetActivityMatrix;
use crate::stats::robust::{extract_geneset_slice, robust_z};

const REQUIRED: &[&str] = &["SRSF_SR", "HNRNP", "U2AF_AXIS"];

pub fn run_stage9(
    activity: &GenesetActivityMatrix,
) -> Result<ExonIntronDefinitionMetrics, InputError> {
    compute(activity)
}

pub fn compute(
    activity: &GenesetActivityMatrix,
) -> Result<ExonIntronDefinitionMetrics, InputError> {
    debug!(genesets = ?activity.genesets, "exon/intron genesets resolved");

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
        return Err(InputError::InsufficientGenesetsExonIntron);
    }

    let start = Instant::now();

    let z_srsf = z_scores.remove("SRSF_SR").unwrap_or_else(|| vec![f32::NAN; n_cells]);
    let z_hnrnp = z_scores.remove("HNRNP").unwrap_or_else(|| vec![f32::NAN; n_cells]);
    let z_u2af = z_scores.remove("U2AF_AXIS").unwrap_or_else(|| vec![f32::NAN; n_cells]);

    let exon_definition_bias: Vec<f32> = (0..n_cells)
        .into_par_iter()
        .map(|cell| {
            let srsf = z_srsf[cell];
            let u2af = z_u2af[cell];
            let hnrnp = z_hnrnp[cell];
            if srsf.is_finite() && u2af.is_finite() && hnrnp.is_finite() {
                (srsf + u2af) - hnrnp
            } else {
                f32::NAN
            }
        })
        .collect();

    info!(
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
