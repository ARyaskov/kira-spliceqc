//! Stage 17: Tier A intron retention index (requires input level L1).

use std::time::Instant;

use tracing::info;

use crate::expression::{ExpressionMatrix, SplicedUnspliced};
use crate::metrics::intron_retention::ExternalIntronRetentionNorms;
use crate::model::intron_retention::IntronRetentionMetrics;
use crate::reference::Strata;
use crate::reference::external::ReferenceFile;

pub fn run_stage17(
    matrix: &dyn ExpressionMatrix,
    layers: &SplicedUnspliced,
    strata: &Strata,
    external: Option<&ReferenceFile>,
) -> IntronRetentionMetrics {
    let start = Instant::now();
    let external_norms = external.map(|file| ExternalIntronRetentionNorms {
        gene_ratios: file.gene_ratios_for(matrix),
        norms: file.intron_retention_norms(),
    });
    let metrics =
        crate::metrics::intron_retention::compute(layers, strata, external_norms.as_ref());
    info!(
        elapsed_ms = start.elapsed().as_millis(),
        undefined_cells = metrics.undefined_cells,
        genes_with_reference = metrics.genes_with_reference,
        high_flags = metrics.intron_retention_high.iter().filter(|f| **f).count(),
        norms = metrics.norm_source,
        "intron retention index computed"
    );
    metrics
}
