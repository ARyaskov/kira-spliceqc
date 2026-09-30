//! Stage 18: cell-cycle phase scores (confounder guard for expression signatures).

use std::time::Instant;

use tracing::info;

use crate::expression::ExpressionMatrix;
use crate::genesets::controls::ControlPool;
use crate::model::cell_cycle::CellCycleMetrics;

pub fn run_stage18(matrix: &dyn ExpressionMatrix, controls: Option<&ControlPool>) -> CellCycleMetrics {
    let start = Instant::now();
    let metrics = crate::metrics::cell_cycle::compute(matrix, controls);
    info!(
        elapsed_ms = start.elapsed().as_millis(),
        "cell-cycle stage complete"
    );
    metrics
}
