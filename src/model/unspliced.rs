/// Tier A per-cell metrics computed from spliced/unspliced count layers.
///
/// All vectors are indexed by cell. Fractions are `NaN` for cells whose
/// spliced + unspliced total is below `min_layer_umis`.
#[derive(Debug, Clone)]
pub struct UnsplicedMetrics {
    /// Provenance of the layers (`mtx-dir:...` or `h5ad-layers:...`).
    pub source: String,
    /// Minimum spliced + unspliced UMIs for a defined fraction.
    pub min_layer_umis: u64,
    pub has_ambiguous: bool,
    pub spliced_umis: Vec<u64>,
    pub unspliced_umis: Vec<u64>,
    pub ambiguous_umis: Vec<u64>,
    /// `U / (S + U)`; ambiguous UMIs are excluded from both terms.
    pub unspliced_fraction: Vec<f32>,
    /// Wilson 95 % interval of `unspliced_fraction`.
    pub unspliced_fraction_ci_low: Vec<f32>,
    pub unspliced_fraction_ci_high: Vec<f32>,
    /// Main-matrix cells that had no column in the layer source.
    pub cells_without_layers: usize,
    /// Cells with an undefined fraction (below `min_layer_umis`).
    pub undefined_cells: usize,
}
