//! Stage 19: Tier B junction metrics (requires input level L2).

use std::time::Instant;

use tracing::info;

use crate::expression::JunctionSet;
use crate::model::junctions::JunctionMetrics;
use crate::reference::Strata;

pub fn run_stage19(junctions: &JunctionSet, strata: &Strata) -> JunctionMetrics {
    let start = Instant::now();
    let metrics = crate::metrics::junctions::compute(junctions, strata);
    info!(
        elapsed_ms = start.elapsed().as_millis(),
        junctions = metrics.n_junctions,
        annotated = metrics.n_annotated,
        cryptic_acceptor_junctions = metrics.n_cryptic_acceptor_junctions,
        skip_junctions = metrics.n_skip_junctions,
        site_groups = metrics.n_site_groups,
        undefined_cells = metrics.undefined_cells,
        cryptic_high = metrics.cryptic_3ss_high.iter().filter(|f| **f).count(),
        skip_high = metrics.exon_skip_high.iter().filter(|f| **f).count(),
        "junction metrics computed"
    );
    metrics
}
