pub trait ExpressionMatrix: Sync {
    fn n_genes(&self) -> usize;
    fn n_cells(&self) -> usize;

    fn gene_symbol(&self, gene_id: usize) -> &str;
    fn cell_name(&self, cell_id: usize) -> &str;

    fn libsize(&self, cell_id: usize) -> u64;
    fn count(&self, gene_id: usize, cell_id: usize) -> u32;

    /// Number of non-zero genes in `cell_id`. Default scans all genes via
    /// `count`; production backends override with O(1).
    fn nnz_cell(&self, cell_id: usize) -> u64 {
        (0..self.n_genes())
            .filter(|&g| self.count(g, cell_id) > 0)
            .count() as u64
    }

    #[inline]
    fn cp10k(&self, gene_id: usize, cell_id: usize) -> f32 {
        let count = self.count(gene_id, cell_id) as f32;
        let libsize = self.libsize(cell_id).max(1) as f32;
        1e4_f32 * count / libsize
    }

    #[inline]
    fn log_cp10k(&self, gene_id: usize, cell_id: usize) -> f32 {
        (1.0 + self.cp10k(gene_id, cell_id)).ln()
    }

    /// Append `log_cp10k(panel[i], cell)` for each `i` in `panel_sorted` order
    /// to `out`. Panel indices must be sorted ascending; sparse backends use a
    /// one-pass merge against the cell's column.
    ///
    /// Default impl performs one `log_cp10k` call per panel gene.
    fn gather_panel_log_cp10k(
        &self,
        panel_sorted: &[u32],
        cell: usize,
        out: &mut Vec<f32>,
    ) {
        out.reserve(panel_sorted.len());
        for &g in panel_sorted {
            out.push(self.log_cp10k(g as usize, cell));
        }
    }

    /// Σ ln(1 + scale * count[g, cell]) over `panel_sorted` (genes absent
    /// from the cell column contribute 0). Used by stage 2 aggregation.
    fn panel_ln1p_scaled_sum(&self, panel_sorted: &[u32], cell: usize, scale: f32) -> f32 {
        let mut sum = 0.0_f32;
        for &g in panel_sorted {
            let c = self.count(g as usize, cell) as f32;
            if c > 0.0 {
                sum += (1.0 + scale * c).ln();
            }
        }
        sum
    }

    /// `(sum, detected)` over `panel_sorted` — used by pipeline_contract panels report.
    fn panel_count_sum_and_detected(
        &self,
        panel_sorted: &[u32],
        cell: usize,
    ) -> (u64, usize) {
        let mut sum = 0u64;
        let mut detected = 0usize;
        for &g in panel_sorted {
            let c = self.count(g as usize, cell);
            if c > 0 {
                detected += 1;
                sum += c as u64;
            }
        }
        (sum, detected)
    }
}

pub mod cache_writer;
pub mod index;
pub mod mmap;

pub use mmap::MmapExpressionMatrix;
