/// Tier B per-cell junction metrics (input level L2).
#[derive(Debug, Clone)]
pub struct JunctionMetrics {
    pub source: String,
    pub n_junctions: usize,
    pub n_annotated: usize,
    /// Unannotated junctions classified as cryptic 3' splice sites.
    pub n_cryptic_acceptor_junctions: usize,
    /// Junctions (annotated or not) that skip at least one annotated exon.
    pub n_skip_junctions: usize,
    /// Alternative splice-site groups (donors with >= 2 acceptors or
    /// acceptors with >= 2 donors) used by `splice_site_shift`.
    pub n_site_groups: usize,
    pub min_junction_umis: u64,
    pub min_ratio_umis: u64,
    pub cells_without_junctions: usize,
    /// Cells below `min_junction_umis` (every ratio undefined).
    pub undefined_cells: usize,

    pub junction_umis: Vec<u64>,
    pub annotated_umis: Vec<u64>,
    /// Unannotated / all junction UMIs.
    pub unannotated_junction_fraction: Vec<f32>,
    /// UMIs on cryptic 3' junctions / (cryptic + canonical partner UMIs).
    pub cryptic_3ss_umis: Vec<u64>,
    pub cryptic_3ss_fraction: Vec<f32>,
    pub cryptic_3ss_fraction_dev: Vec<f32>,
    pub cryptic_3ss_high: Vec<bool>,
    /// Skip-junction UMIs / (skip + inclusion UMIs).
    pub exon_skip_umis: Vec<u64>,
    pub exon_skip_fraction: Vec<f32>,
    pub exon_skip_fraction_dev: Vec<f32>,
    pub exon_skip_high: Vec<bool>,
    /// SpliZ-like usage shift: median over alternative-site groups of the
    /// standardized rank deviation (|z| / 0.6745; ~1 under the null).
    pub splice_site_shift: Vec<f32>,
    pub splice_site_shift_dev: Vec<f32>,
    pub splice_site_shift_high: Vec<bool>,
    pub site_groups_used: Vec<u32>,

    pub cryptic_reference: Vec<crate::reference::StratumStat>,
    pub skip_reference: Vec<crate::reference::StratumStat>,
    pub shift_reference: Vec<crate::reference::StratumStat>,
}
