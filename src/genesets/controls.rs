//! Control-gene background for expression panel scores.
//!
//! A panel mean of `log1p(cp10k)` over a handful of genes rises with library
//! size and detection rate on sparse data, so it mostly measures depth. The
//! standard remedy (Tirosh et al. 2016 Science; Seurat `AddModuleScore`) is
//! to subtract the mean of control genes with the same expression level:
//! for every panel gene the pool takes the `CONTROLS_PER_GENE` eligible genes
//! nearest to it in mean expression, and the panel's control set is their
//! union. Genes that belong to any splicing panel (or to a caller-supplied
//! exclusion list, e.g. cell-cycle genes) are never controls.

use std::collections::BTreeSet;

/// Controls drawn per panel gene (Seurat default: 100 per bin of 24; 50
/// nearest neighbours gives a comparable pool for 3-40-gene panels).
pub const CONTROLS_PER_GENE: usize = 50;

#[derive(Debug, Clone)]
pub struct ControlPool {
    /// Eligible gene indices sorted by mean expression (ascending).
    ranked: Vec<u32>,
    /// Mean expression per gene (all genes, indexed by gene id).
    means: Vec<f32>,
    pub controls_per_gene: usize,
    /// Genes excluded from the pool (panel members, caller exclusions, unexpressed).
    pub excluded: usize,
}

impl ControlPool {
    /// `means` = per-gene mean of `log1p(cp10k)` over cells; `exclude` = gene
    /// ids that must not serve as controls (any order, duplicates allowed).
    pub fn new(means: Vec<f32>, exclude: &[u32]) -> Self {
        let excluded: BTreeSet<u32> = exclude.iter().copied().collect();
        let mut ranked: Vec<u32> = (0..means.len() as u32)
            .filter(|g| !excluded.contains(g) && means[*g as usize] > 0.0 && means[*g as usize].is_finite())
            .collect();
        ranked.sort_by(|a, b| {
            means[*a as usize]
                .partial_cmp(&means[*b as usize])
                .unwrap()
                .then_with(|| a.cmp(b))
        });
        let excluded = means.len() - ranked.len();
        Self {
            ranked,
            means,
            controls_per_gene: CONTROLS_PER_GENE,
            excluded,
        }
    }

    pub fn n_eligible(&self) -> usize {
        self.ranked.len()
    }

    /// Sorted, de-duplicated control set for a panel: the union over panel
    /// genes of the `controls_per_gene` eligible genes nearest in mean
    /// expression. Empty when the pool has no eligible genes.
    pub fn controls_for(&self, panel: &[u32]) -> Vec<u32> {
        if self.ranked.is_empty() {
            return Vec::new();
        }
        let panel_set: BTreeSet<u32> = panel.iter().copied().collect();
        let mut set: BTreeSet<u32> = BTreeSet::new();
        for &g in panel {
            let target = self.means[g as usize];
            // Insertion point of the panel gene's mean in the ranked pool.
            let pos = self
                .ranked
                .partition_point(|&r| self.means[r as usize] < target);
            // Expand a window around `pos`, always taking the closer side;
            // the panel's own genes never count as controls.
            let (mut lo, mut hi) = (pos, pos);
            let mut taken = 0usize;
            while taken < self.controls_per_gene && (lo > 0 || hi < self.ranked.len()) {
                let take_left = if lo == 0 {
                    false
                } else if hi >= self.ranked.len() {
                    true
                } else {
                    (target - self.means[self.ranked[lo - 1] as usize]).abs()
                        <= (self.means[self.ranked[hi] as usize] - target).abs()
                };
                let candidate = if take_left {
                    lo -= 1;
                    self.ranked[lo]
                } else {
                    hi += 1;
                    self.ranked[hi - 1]
                };
                if !panel_set.contains(&candidate) {
                    set.insert(candidate);
                    taken += 1;
                }
            }
        }
        set.into_iter().collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nearest_by_mean_excluding_panel_and_listed_genes() {
        // Means 0.0 (unexpressed), then 0.1 .. 2.0 in steps of 0.1 for genes 1..=20.
        let mut means = vec![0.0f32];
        means.extend((1..=20).map(|i| i as f32 * 0.1));
        let pool = ControlPool::new(means.clone(), &[10, 11]);
        assert_eq!(pool.n_eligible(), 18); // 20 expressed minus 2 excluded
        assert_eq!(pool.excluded, 3);

        let mut small = pool.clone();
        small.controls_per_gene = 4;
        // Panel gene 10 (mean 1.0): nearest eligible are 9, 12, 8, 13 (ties resolved toward the closer side).
        let ctrl = small.controls_for(&[10]);
        assert_eq!(ctrl.len(), 4);
        assert!(!ctrl.contains(&10) && !ctrl.contains(&11) && !ctrl.contains(&0));
        assert!(ctrl.contains(&9) && ctrl.contains(&12));
        // A panel gene that is still in the pool is never its own control.
        let pool_with_panel_gene = ControlPool::new(means.clone(), &[]);
        let mut p = pool_with_panel_gene;
        p.controls_per_gene = 2;
        let ctrl = p.controls_for(&[10]);
        assert_eq!(ctrl.len(), 2);
        assert!(!ctrl.contains(&10));
        // Union over two panel genes is de-duplicated and sorted.
        let ctrl = small.controls_for(&[10, 11]);
        assert!(ctrl.windows(2).all(|w| w[0] < w[1]));
        assert!(ctrl.len() <= 8);
    }

    #[test]
    fn window_clamps_at_pool_edges() {
        let means: Vec<f32> = (0..10).map(|i| 0.1 + i as f32 * 0.1).collect();
        let mut pool = ControlPool::new(means, &[]);
        pool.controls_per_gene = 3;
        // Lowest gene: window extends upward only.
        assert_eq!(pool.controls_for(&[0]), vec![1, 2, 3]);
        // Highest gene: window extends downward only.
        assert_eq!(pool.controls_for(&[9]), vec![6, 7, 8]);
        // Pool smaller than the request: everything eligible.
        pool.controls_per_gene = 100;
        assert_eq!(pool.controls_for(&[5]).len(), 9);
    }

    #[test]
    fn empty_pool_gives_no_controls() {
        let pool = ControlPool::new(vec![0.0; 5], &[]);
        assert!(pool.controls_for(&[1, 2]).is_empty());
    }
}
