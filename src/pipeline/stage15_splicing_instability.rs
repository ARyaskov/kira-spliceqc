use crate::expression::ExpressionMatrix;
use crate::expression::MmapExpressionMatrix;
use crate::model::splicing_instability::SplicingInstabilityMetrics;

pub fn run_stage15(matrix: &MmapExpressionMatrix) -> SplicingInstabilityMetrics {
    compute(matrix)
}

pub fn compute(matrix: &dyn ExpressionMatrix) -> SplicingInstabilityMetrics {
    crate::metrics::splicing_instability::compute(matrix)
}
