//! Stage 17: Tier A intron retention index (requires input level L1).

use std::time::Instant;

use tracing::info;

use crate::expression::SplicedUnspliced;
use crate::model::intron_retention::IntronRetentionMetrics;
use crate::reference::Strata;

pub fn run_stage17(layers: &SplicedUnspliced, strata: &Strata) -> IntronRetentionMetrics {
    let start = Instant::now();
    let metrics = crate::metrics::intron_retention::compute(layers, strata);
    info!(
        elapsed_ms = start.elapsed().as_millis(),
        undefined_cells = metrics.undefined_cells,
        genes_with_reference = metrics.genes_with_reference,
        high_flags = metrics.intron_retention_high.iter().filter(|f| **f).count(),
        "intron retention index computed"
    );
    metrics
}
