#[derive(Debug, Clone)]
pub struct SplicingInstabilityMetrics {
    pub panel_version: &'static str,
    pub min_genes: usize,
    pub conflict_panel_enabled: bool,
    pub nmd_panel_enabled: bool,
    pub splice_core: Vec<f32>,
    pub rbp_core: Vec<f32>,
    pub rloop_resolve_core: Vec<f32>,
    pub conflict_risk_core: Vec<f32>,
    pub nmd_core: Vec<f32>,
    pub sos: Vec<f32>,
    pub rlr: Vec<f32>,
    pub sii: Vec<f32>,
    /// Stratum-standardized signed composites behind the flags.
    pub sos_dev: Vec<f32>,
    pub rlr_dev: Vec<f32>,
    pub sii_dev: Vec<f32>,
    pub splice_overload_high: Vec<bool>,
    pub rloop_risk_high: Vec<bool>,
    pub splicing_instability_high: Vec<bool>,
    pub genome_instability_splicing_flag: Vec<bool>,
    pub z_reference: SplicingInstabilityZReference,
    pub global_stats: SplicingInstabilityGlobalStats,
    pub cluster_stats: Vec<SplicingInstabilityClusterStat>,
    pub missingness: SplicingInstabilityMissingness,
}

#[derive(Debug, Clone)]
pub struct SplicingInstabilityRobustRef {
    pub median: f32,
    pub mad: f32,
}

#[derive(Debug, Clone)]
pub struct SplicingInstabilityZReference {
    pub splice_core: SplicingInstabilityRobustRef,
    pub rbp_core: SplicingInstabilityRobustRef,
    pub rloop_resolve_core: SplicingInstabilityRobustRef,
    pub conflict_risk_core: Option<SplicingInstabilityRobustRef>,
    pub nmd_core: Option<SplicingInstabilityRobustRef>,
}

#[derive(Debug, Clone)]
pub struct SplicingInstabilityGlobalStats {
    pub sos_p50: f32,
    pub sos_p90: f32,
    pub rlr_p50: f32,
    pub rlr_p90: f32,
    pub sii_p50: f32,
    pub sii_p90: f32,
}

#[derive(Debug, Clone)]
pub struct SplicingInstabilityClusterStat {
    pub cluster_id: String,
    pub n_cells: usize,
    pub sos_median: f32,
    pub sos_p10: f32,
    pub sos_p90: f32,
    pub rlr_median: f32,
    pub rlr_p10: f32,
    pub rlr_p90: f32,
    pub sii_median: f32,
    pub sii_p10: f32,
    pub sii_p90: f32,
    pub splice_overload_high_fraction: f32,
    pub rloop_risk_high_fraction: f32,
    pub splicing_instability_high_fraction: f32,
    pub genome_instability_splicing_flag_fraction: f32,
}

#[derive(Debug, Clone)]
pub struct PanelCoverage {
    pub panel_name: &'static str,
    pub genes_defined: usize,
    pub genes_mapped: usize,
}

#[derive(Debug, Clone)]
pub struct SplicingInstabilityMissingness {
    pub splice_core_nan_cells: usize,
    pub rbp_core_nan_cells: usize,
    pub rloop_resolve_core_nan_cells: usize,
    pub conflict_risk_core_nan_cells: usize,
    pub nmd_core_nan_cells: usize,
    pub sos_nan_cells: usize,
    pub rlr_nan_cells: usize,
    pub sii_nan_cells: usize,
    pub panel_coverage: Vec<PanelCoverage>,
}

/// Flag rule shared with every other outlier flag (see `reference`): the
/// signed, stratum-standardized composite must deviate by at least
/// `DEVIATION_THRESHOLD` with a BH-adjusted p below `FLAG_FDR`. The fixed
/// cut-offs of v0.3 (SOS >= 2.0, RLR >= 1.5, SII >= 2.0) fired on 4 % of
/// null cells and are gone.
pub const COMPOSITE_FLAG_RULE: &str =
    "signed composite dev >= 3 and BH-adjusted p < 0.05 within stratum";
