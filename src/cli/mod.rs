pub mod config;
pub mod run;

use crate::expression::MmapExpressionMatrix;
use crate::input::InputDescriptor;
use crate::model::assembly_phase::AssemblyPhaseImbalanceMetrics;
use crate::model::collapse::SpliceosomeCollapseMetrics;
use crate::model::coupling::CouplingStressMetrics;
use crate::model::cryptic_risk::CrypticSplicingRiskMetrics;
use crate::model::exon_intron_bias::ExonIntronDefinitionMetrics;
use crate::model::geneset_activity::GenesetActivityMatrix;
use crate::model::imbalance::SpliceosomeImbalanceMetrics;
use crate::model::isoform_dispersion::IsoformDispersionMetrics;
use crate::model::missplicing::MissplicingMetrics;
use crate::model::sis::SpliceIntegrityMetrics;
use crate::model::splicing_noise::SplicingNoiseMetrics;
use crate::model::timecourse::TimecourseSplicingMetrics;

pub struct PipelineContext {
    pub stage0: InputDescriptor,
    pub stage1: MmapExpressionMatrix,
    pub stage2: GenesetActivityMatrix,
    pub stage3: IsoformDispersionMetrics,
    pub stage4: MissplicingMetrics,
    pub stage5: SpliceosomeImbalanceMetrics,
    pub stage6: SpliceIntegrityMetrics,
    pub stage8: Option<CouplingStressMetrics>,
    pub stage9: Option<ExonIntronDefinitionMetrics>,
    pub stage10: Option<AssemblyPhaseImbalanceMetrics>,
    pub stage11: Option<SplicingNoiseMetrics>,
    pub stage12: Option<CrypticSplicingRiskMetrics>,
    pub stage13: Option<SpliceosomeCollapseMetrics>,
    pub stage14: Option<TimecourseSplicingMetrics>,
}
