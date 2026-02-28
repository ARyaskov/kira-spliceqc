pub trait ExpressionMatrix: Sync {
    fn n_genes(&self) -> usize;
    fn n_cells(&self) -> usize;

    fn gene_symbol(&self, gene_id: usize) -> &str;
    fn cell_name(&self, cell_id: usize) -> &str;

    fn libsize(&self, cell_id: usize) -> u64;
    fn count(&self, gene_id: usize, cell_id: usize) -> u32;
    fn nnz_cell(&self, cell_id: usize) -> u64 {
        let mut nnz = 0u64;
        for gene_id in 0..self.n_genes() {
            if self.count(gene_id, cell_id) > 0 {
                nnz += 1;
            }
        }
        nnz
    }

    fn cp10k(&self, gene_id: usize, cell_id: usize) -> f32 {
        let count = self.count(gene_id, cell_id) as f32;
        let libsize = self.libsize(cell_id).max(1) as f32;
        1e4_f32 * count / libsize
    }

    fn log_cp10k(&self, gene_id: usize, cell_id: usize) -> f32 {
        (1.0 + self.cp10k(gene_id, cell_id)).ln()
    }
}

pub mod cache_writer;
pub mod index;
pub mod mmap;

pub use mmap::MmapExpressionMatrix;
