//! Spliced / unspliced / ambiguous count layers (input level L1).
//!
//! Layers are stored in memory as CSC (one column per cell, rows sorted by
//! gene) in the *same* gene and cell order as the main expression matrix, so
//! per-cell sparse merges between the L0 matrix and any layer are direct.

/// One count layer, CSC with rows sorted ascending within every column.
#[derive(Debug, Clone)]
pub struct LayerMatrix {
    n_genes: usize,
    n_cells: usize,
    col_ptr: Vec<u64>,
    row_idx: Vec<u32>,
    values: Vec<u32>,
}

impl LayerMatrix {
    /// Builds a layer from `(gene, cell, count)` triplets. Triplets may be in
    /// any order; duplicates of the same `(gene, cell)` are summed.
    pub fn from_triplets(
        n_genes: usize,
        n_cells: usize,
        mut triplets: Vec<(u32, u32, u32)>,
    ) -> Self {
        triplets.sort_unstable_by_key(|&(gene, cell, _)| (cell, gene));

        let mut col_ptr: Vec<u64> = Vec::with_capacity(n_cells + 1);
        let mut row_idx: Vec<u32> = Vec::with_capacity(triplets.len());
        let mut values: Vec<u32> = Vec::with_capacity(triplets.len());
        col_ptr.push(0);
        let mut cur_col: u32 = 0;
        for (gene, cell, count) in triplets {
            if count == 0 {
                continue;
            }
            while cur_col < cell {
                col_ptr.push(row_idx.len() as u64);
                cur_col += 1;
            }
            // Same (gene, cell) as the previous entry of *this* column: sum.
            let col_start = *col_ptr.last().unwrap() as usize;
            if row_idx.len() > col_start && row_idx[row_idx.len() - 1] == gene {
                let last = values.len() - 1;
                values[last] = values[last].saturating_add(count);
                continue;
            }
            row_idx.push(gene);
            values.push(count);
        }
        while col_ptr.len() <= n_cells {
            col_ptr.push(row_idx.len() as u64);
        }

        Self {
            n_genes,
            n_cells,
            col_ptr,
            row_idx,
            values,
        }
    }

    pub fn n_genes(&self) -> usize {
        self.n_genes
    }

    pub fn n_cells(&self) -> usize {
        self.n_cells
    }

    pub fn nnz(&self) -> usize {
        self.values.len()
    }

    /// `(gene indices, counts)` of one cell's column.
    #[inline]
    pub fn cell(&self, cell: usize) -> (&[u32], &[u32]) {
        let s = self.col_ptr[cell] as usize;
        let e = self.col_ptr[cell + 1] as usize;
        (&self.row_idx[s..e], &self.values[s..e])
    }

    /// Sum of counts in one cell.
    pub fn cell_total(&self, cell: usize) -> u64 {
        self.cell(cell).1.iter().map(|&v| v as u64).sum()
    }

    /// Count for one `(gene, cell)`; 0 when absent.
    pub fn count(&self, gene: usize, cell: usize) -> u32 {
        let (rows, vals) = self.cell(cell);
        match rows.binary_search(&(gene as u32)) {
            Ok(i) => vals[i],
            Err(_) => 0,
        }
    }
}

/// The spliced / unspliced (/ ambiguous) layer set of one dataset.
#[derive(Debug, Clone)]
pub struct SplicedUnspliced {
    pub spliced: LayerMatrix,
    pub unspliced: LayerMatrix,
    pub ambiguous: Option<LayerMatrix>,
    /// Where the layers came from (for provenance).
    pub source: String,
    /// Cells of the main matrix that had no column in the layer source
    /// (their layer counts are zero).
    pub cells_without_layers: usize,
}

impl SplicedUnspliced {
    pub fn n_cells(&self) -> usize {
        self.spliced.n_cells()
    }

    pub fn n_genes(&self) -> usize {
        self.spliced.n_genes()
    }

    /// Visits every gene that has a spliced or unspliced count in `cell`,
    /// as `(gene, spliced, unspliced)`; genes are visited in ascending order.
    pub fn for_each_gene(&self, cell: usize, mut f: impl FnMut(u32, u32, u32)) {
        let (sg, sv) = self.spliced.cell(cell);
        let (ug, uv) = self.unspliced.cell(cell);
        let (mut i, mut j) = (0usize, 0usize);
        while i < sg.len() || j < ug.len() {
            match (sg.get(i), ug.get(j)) {
                (Some(&g1), Some(&g2)) if g1 == g2 => {
                    f(g1, sv[i], uv[j]);
                    i += 1;
                    j += 1;
                }
                (Some(&g1), Some(&g2)) if g1 < g2 => {
                    f(g1, sv[i], 0);
                    i += 1;
                }
                (Some(_), Some(&g2)) => {
                    f(g2, 0, uv[j]);
                    j += 1;
                }
                (Some(&g1), None) => {
                    f(g1, sv[i], 0);
                    i += 1;
                }
                (None, Some(&g2)) => {
                    f(g2, 0, uv[j]);
                    j += 1;
                }
                (None, None) => unreachable!(),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_triplets_sorts_and_merges_duplicates() {
        let m = LayerMatrix::from_triplets(3, 2, vec![(2, 1, 4), (0, 0, 1), (2, 1, 1), (1, 0, 2)]);
        assert_eq!(m.cell(0), (&[0u32, 1][..], &[1u32, 2][..]));
        assert_eq!(m.cell(1), (&[2u32][..], &[5u32][..]));
        assert_eq!(m.cell_total(1), 5);
        assert_eq!(m.count(1, 0), 2);
        assert_eq!(m.count(1, 1), 0);
        assert_eq!(m.nnz(), 3);
    }

    #[test]
    fn same_gene_in_adjacent_columns_is_not_merged() {
        let m = LayerMatrix::from_triplets(3, 2, vec![(2, 0, 1), (2, 1, 4)]);
        assert_eq!(m.cell(0), (&[2u32][..], &[1u32][..]));
        assert_eq!(m.cell(1), (&[2u32][..], &[4u32][..]));
    }

    #[test]
    fn empty_trailing_cells_get_empty_columns() {
        let m = LayerMatrix::from_triplets(2, 4, vec![(0, 0, 3)]);
        assert_eq!(m.cell(3).0.len(), 0);
        assert_eq!(m.cell_total(0), 3);
    }

    #[test]
    fn for_each_gene_merges_layers() {
        let su = SplicedUnspliced {
            spliced: LayerMatrix::from_triplets(4, 1, vec![(0, 0, 5), (2, 0, 1)]),
            unspliced: LayerMatrix::from_triplets(4, 1, vec![(2, 0, 2), (3, 0, 7)]),
            ambiguous: None,
            source: "test".to_string(),
            cells_without_layers: 0,
        };
        let mut seen = Vec::new();
        su.for_each_gene(0, |g, s, u| seen.push((g, s, u)));
        assert_eq!(seen, vec![(0, 5, 0), (2, 1, 2), (3, 0, 7)]);
    }
}
