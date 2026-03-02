//! Optional junction-aware instability metrics are not yet wired in runtime input modes.
//! TODO: enable JE/IRB/CSP once per-cell junction and intron-count tables are available.

#[derive(Debug, Clone)]
pub struct JunctionAwareInstabilityMetrics {
    pub junction_entropy: Vec<f32>,
    pub intron_retention_burden: Vec<f32>,
    pub cryptic_splicing_proxy: Vec<f32>,
}
