//! Junction count matrix (input level L2): one row per splice junction, one
//! column per cell, in the main matrix's cell order.

use crate::expression::LayerMatrix;
use crate::io::junctions::Junction;

#[derive(Debug, Clone)]
pub struct JunctionSet {
    pub junctions: Vec<Junction>,
    /// `n_junctions x n_cells` CSC counts (rows = junction ids).
    pub counts: LayerMatrix,
    pub source: String,
    pub cells_without_junctions: usize,
}

impl JunctionSet {
    pub fn n_junctions(&self) -> usize {
        self.junctions.len()
    }

    pub fn n_cells(&self) -> usize {
        self.counts.n_cells()
    }
}
