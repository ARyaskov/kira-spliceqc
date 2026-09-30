//! Stage 16: Tier A unspliced fraction (requires input level L1).

use std::time::Instant;

use tracing::info;

use crate::expression::SplicedUnspliced;
use crate::model::unspliced::UnsplicedMetrics;

pub fn run_stage16(layers: &SplicedUnspliced) -> UnsplicedMetrics {
    let start = Instant::now();
    let metrics = crate::metrics::unspliced::compute(layers);
    info!(
        elapsed_ms = start.elapsed().as_millis(),
        undefined_cells = metrics.undefined_cells,
        cells_without_layers = metrics.cells_without_layers,
        "unspliced fraction computed"
    );
    metrics
}
