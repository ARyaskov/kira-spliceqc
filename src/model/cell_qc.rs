/// Per-cell QC flags that mark cells whose metrics should not be
/// interpreted and that are excluded from reference-norm computation.
#[derive(Debug, Clone)]
pub struct CellQc {
    pub min_counts: u64,
    pub min_genes: u64,
    /// `libsize < min_counts` or `nnz < min_genes`.
    pub low_depth: Vec<bool>,
    /// Marked as a doublet by a metadata column (`predicted_doublet`, ...).
    pub doublet: Vec<bool>,
    /// Metadata column the doublet flag came from, if any.
    pub doublet_column: Option<String>,
}

impl CellQc {
    pub fn none(n_cells: usize) -> Self {
        Self {
            min_counts: 0,
            min_genes: 0,
            low_depth: vec![false; n_cells],
            doublet: vec![false; n_cells],
            doublet_column: None,
        }
    }

    pub fn n_cells(&self) -> usize {
        self.low_depth.len()
    }

    /// Cells excluded from reference norms (low depth or doublet).
    pub fn excluded(&self) -> Vec<bool> {
        self.low_depth
            .iter()
            .zip(&self.doublet)
            .map(|(l, d)| *l || *d)
            .collect()
    }

    pub fn n_low_depth(&self) -> usize {
        self.low_depth.iter().filter(|f| **f).count()
    }

    pub fn n_doublet(&self) -> usize {
        self.doublet.iter().filter(|f| **f).count()
    }
}
