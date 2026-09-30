/// Seurat-style cell-cycle phase call.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CellCyclePhase {
    G1,
    S,
    G2M,
    /// Scores undefined (too few genes mapped).
    Unknown,
}

impl CellCyclePhase {
    pub fn as_str(&self) -> &'static str {
        match self {
            CellCyclePhase::G1 => "G1",
            CellCyclePhase::S => "S",
            CellCyclePhase::G2M => "G2M",
            CellCyclePhase::Unknown => "",
        }
    }
}

/// Per-cell cell-cycle scores (confounder annotation).
#[derive(Debug, Clone)]
pub struct CellCycleMetrics {
    pub s_genes_mapped: usize,
    pub g2m_genes_mapped: usize,
    /// Control-corrected mean log1p(cp10k) of the S-phase genes.
    pub s_score: Vec<f32>,
    /// Control-corrected mean log1p(cp10k) of the G2/M genes.
    pub g2m_score: Vec<f32>,
    pub phase: Vec<CellCyclePhase>,
    /// Phase is S or G2M.
    pub cycling: Vec<bool>,
}
