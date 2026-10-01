use crate::expression::ExpressionMatrix;
use crate::expression::MmapExpressionMatrix;
use crate::genesets::controls::ControlPool;
use crate::model::splicing_instability::SplicingInstabilityMetrics;
use crate::reference::Strata;
use crate::reference::external::ReferenceFile;

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

pub fn compute_with(
    matrix: &dyn ExpressionMatrix,
    controls: Option<&ControlPool>,
    strata: &Strata,
    external: Option<&ReferenceFile>,
) -> SplicingInstabilityMetrics {
    crate::metrics::splicing_instability::compute_with(matrix, controls, strata, external)
}
