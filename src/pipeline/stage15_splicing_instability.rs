use crate::expression::ExpressionMatrix;
use crate::expression::MmapExpressionMatrix;
use crate::genesets::controls::ControlPool;
use crate::model::splicing_instability::SplicingInstabilityMetrics;
use crate::reference::Strata;

pub fn run_stage15(matrix: &MmapExpressionMatrix) -> SplicingInstabilityMetrics {
    compute(matrix, None, &Strata::global(matrix.n_cells()))
}

pub fn compute(
    matrix: &dyn ExpressionMatrix,
    controls: Option<&ControlPool>,
    strata: &Strata,
) -> SplicingInstabilityMetrics {
    crate::metrics::splicing_instability::compute(matrix, controls, strata)
}
