/// Tier A per-cell intron retention index from spliced/unspliced layers.
#[derive(Debug, Clone)]
pub struct IntronRetentionMetrics {
    /// Minimum `S + U` UMIs of a gene in a cell for the gene to count.
    pub min_gene_umis: u32,
    /// Minimum genes with a defined ratio for a defined index.
    pub min_genes: usize,
    /// Median over genes of `log2(IR_gc / IR_g,ref)`; NaN when undefined.
    pub intron_retention_index: Vec<f32>,
    /// Robust z-score of the index within the cell's reference stratum.
    pub intron_retention_index_dev: Vec<f32>,
    /// MAD over genes of the per-gene log2 ratios (gene-specific vs global shift).
    pub ir_gene_dispersion: Vec<f32>,
    /// Genes that contributed to the index.
    pub ir_genes_used: Vec<u32>,
    /// Index far above its stratum (`dev >= 3`, BH-adjusted p < 0.05).
    pub intron_retention_high: Vec<bool>,
    /// Cells with an undefined index.
    pub undefined_cells: usize,
    /// Per-stratum reference of the index (median / MAD).
    pub reference: Vec<crate::reference::StratumStat>,
    /// Genes with a defined reference ratio in at least one stratum.
    pub genes_with_reference: usize,
}
