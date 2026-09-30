//! Stage 16: Tier A unspliced fraction (requires input level L1).

use std::time::Instant;

use tracing::info;

use crate::expression::SplicedUnspliced;
use crate::model::unspliced::UnsplicedMetrics;
use crate::reference::Strata;
use crate::reference::external::ReferenceFile;

pub fn run_stage16(
    layers: &SplicedUnspliced,
    strata: &Strata,
    external: Option<&ReferenceFile>,
) -> UnsplicedMetrics {
    let start = Instant::now();
    let metrics = crate::metrics::unspliced::compute(layers, strata, external);
    info!(
        elapsed_ms = start.elapsed().as_millis(),
        undefined_cells = metrics.undefined_cells,
        cells_without_layers = metrics.cells_without_layers,
        nuclear_fraction_flags = metrics.nuclear_fraction_flag.iter().filter(|f| **f).count(),
        reference = strata.mode.as_str(),
        norms = metrics.norm_source,
        "unspliced fraction computed"
    );
    metrics
}
