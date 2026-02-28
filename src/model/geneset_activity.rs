#[derive(Debug, Clone)]
pub struct GenesetActivityMatrix {
    pub genesets: Vec<String>,
    pub axes: Vec<String>,
    pub values: Vec<f32>,
    pub n_cells: usize,
}

impl GenesetActivityMatrix {
    pub fn value(&self, geneset_idx: usize, cell_id: usize) -> f32 {
        let index = geneset_idx * self.n_cells + cell_id;
        self.values[index]
    }
}
