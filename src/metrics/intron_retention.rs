//! Tier A: per-cell intron retention index (IRI).
//!
//! For gene `g` in cell `c` the unspliced ratio is shrunk toward its
//! reference stratum with a beta prior (empirical Bayes):
//!
//! `IR_gc = (U_gc + a_gs) / (S_gc + U_gc + a_gs + b_gs)`,
//! `a_gs = K * p_gs`, `b_gs = K * (1 - p_gs)`,
//!
//! where `p_gs` is the pooled unspliced ratio of gene `g` over the cells of
//! stratum `s` and `K = PRIOR_STRENGTH` pseudo-counts. The prior is fixed and
//! modest on purpose: a method-of-moments prior strength fitted to the
//! population collapses to hundreds of pseudo-counts whenever most cells
//! agree, and then shrinks exactly the outlier cells the index is meant to
//! reveal. With `K = 10`, a gene with 5 UMIs is pulled two thirds of the way
//! to the reference while a gene with 50 UMIs keeps five sixths of its own
//! signal. The cell's index is the median over informative genes of
//! `log2(IR_gc / p_gs)`; the MAD of the same values separates a global shift
//! (retention in every gene) from gene-specific changes.
//!
//! Adapted from the intron-retention ratio of IRFinder (Middleton et al. 2017)
//! to sparse per-cell counts.

use rayon::prelude::*;

use crate::expression::SplicedUnspliced;
use crate::model::intron_retention::IntronRetentionMetrics;
use crate::reference::{Strata, flag_outliers, robust_z_by_stratum};
use crate::stats::robust::{mad, median};

/// A gene counts for a cell when it has at least this many `S + U` UMIs.
pub const MIN_GENE_UMIS: u32 = 5;
/// A cell's index is defined when at least this many genes count.
pub const MIN_GENES: usize = 20;
/// A gene has a reference in a stratum when at least this many of the
/// stratum's cells carry >= `MIN_GENE_UMIS` of it.
pub const MIN_CELLS_PER_GENE: usize = 10;
/// Prior strength (pseudo-counts) of the beta shrinkage toward the stratum ratio.
pub const PRIOR_STRENGTH: f64 = 10.0;

/// Per-gene reference of one stratum.
struct GeneReference {
    /// Pooled unspliced ratio `p_gs`; NaN when undefined.
    p: Vec<f64>,
}

pub fn compute(layers: &SplicedUnspliced, strata: &Strata) -> IntronRetentionMetrics {
    let n_cells = layers.n_cells();
    let n_genes = layers.n_genes();
    debug_assert_eq!(strata.n_cells(), n_cells);

    let references: Vec<GeneReference> = strata
        .members()
        .par_iter()
        .map(|members| gene_reference(layers, members, n_genes))
        .collect();
    let genes_with_reference = (0..n_genes)
        .filter(|&g| references.iter().any(|r| r.p[g].is_finite()))
        .count();

    let rows: Vec<(f32, f32, u32)> = (0..n_cells)
        .into_par_iter()
        .map_init(Vec::<f32>::new, |ratios, cell| {
            ratios.clear();
            let reference = &references[strata.labels[cell] as usize];
            layers.for_each_gene(cell, |gene, s, u| {
                let n = s + u;
                if n < MIN_GENE_UMIS {
                    return;
                }
                let g = gene as usize;
                let p = reference.p[g];
                if !p.is_finite() || p <= 0.0 || p >= 1.0 {
                    return;
                }
                let a = PRIOR_STRENGTH * p;
                let b = PRIOR_STRENGTH * (1.0 - p);
                let ir = (u as f64 + a) / (n as f64 + a + b);
                ratios.push((ir / p).log2() as f32);
            });
            let used = ratios.len() as u32;
            if ratios.len() < MIN_GENES {
                (f32::NAN, f32::NAN, used)
            } else {
                let med = median(ratios);
                (med, mad(ratios, med), used)
            }
        })
        .collect();

    let mut index = Vec::with_capacity(n_cells);
    let mut dispersion = Vec::with_capacity(n_cells);
    let mut genes_used = Vec::with_capacity(n_cells);
    let mut undefined_cells = 0usize;
    for (i, d, u) in rows {
        if !i.is_finite() {
            undefined_cells += 1;
        }
        index.push(i);
        dispersion.push(d);
        genes_used.push(u);
    }

    let (dev, reference) = robust_z_by_stratum(&index, strata);
    let high = flag_outliers(&dev, strata, 1.0);

    IntronRetentionMetrics {
        min_gene_umis: MIN_GENE_UMIS,
        min_genes: MIN_GENES,
        intron_retention_index: index,
        intron_retention_index_dev: dev,
        ir_gene_dispersion: dispersion,
        ir_genes_used: genes_used,
        intron_retention_high: high,
        undefined_cells,
        reference,
        genes_with_reference,
    }
}

/// Pooled unspliced ratio per gene over `members` (NaN when fewer than
/// `MIN_CELLS_PER_GENE` cells carry the gene).
fn gene_reference(layers: &SplicedUnspliced, members: &[usize], n_genes: usize) -> GeneReference {
    let mut sum_u = vec![0f64; n_genes];
    let mut sum_n = vec![0f64; n_genes];
    let mut count = vec![0u32; n_genes];
    for &cell in members {
        layers.for_each_gene(cell, |gene, s, u| {
            let n = s + u;
            if n < MIN_GENE_UMIS {
                return;
            }
            let g = gene as usize;
            sum_u[g] += u as f64;
            sum_n[g] += n as f64;
            count[g] += 1;
        });
    }

    let mut p_ref = vec![f64::NAN; n_genes];
    for g in 0..n_genes {
        if (count[g] as usize) < MIN_CELLS_PER_GENE || sum_n[g] <= 0.0 {
            continue;
        }
        // Genes never (or always) unspliced carry no ratio information; a
        // half-count floor keeps them defined but shrunk to the boundary.
        let floor = 0.5 / sum_n[g];
        p_ref[g] = (sum_u[g] / sum_n[g]).clamp(floor, 1.0 - floor);
    }
    GeneReference { p: p_ref }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expression::LayerMatrix;

    /// `n_cells` cells x `n_genes` genes, every gene with `s` spliced and `u`
    /// unspliced UMIs; `boost` cells get `u` doubled.
    fn layers(n_cells: usize, n_genes: usize, s: u32, u: u32, boost: &[usize]) -> SplicedUnspliced {
        let mut sp = Vec::new();
        let mut un = Vec::new();
        for c in 0..n_cells {
            let factor = if boost.contains(&c) { 2 } else { 1 };
            for g in 0..n_genes {
                // small deterministic jitter so per-cell ratios vary
                let jitter = ((c + g) % 3) as u32;
                sp.push((g as u32, c as u32, s + jitter));
                un.push((g as u32, c as u32, u * factor + jitter));
            }
        }
        SplicedUnspliced {
            spliced: LayerMatrix::from_triplets(n_genes, n_cells, sp),
            unspliced: LayerMatrix::from_triplets(n_genes, n_cells, un),
            ambiguous: None,
            source: "test".to_string(),
            cells_without_layers: 0,
        }
    }

    #[test]
    fn boosted_cells_have_positive_index_and_are_flagged() {
        let boosted = [3usize, 40];
        let l = layers(80, 30, 30, 10, &boosted);
        let m = compute(&l, &Strata::global(80));
        assert_eq!(m.undefined_cells, 0);
        assert_eq!(m.genes_with_reference, 30);
        // Typical cell: ratio at the reference -> index ~ 0.
        assert!(m.intron_retention_index[0].abs() < 0.15, "{}", m.intron_retention_index[0]);
        // Boosted cells: unspliced doubled -> index clearly positive (shrinkage
        // pulls it below the naive log2(2) but keeps the sign and order).
        for &c in &boosted {
            assert!(m.intron_retention_index[c] > 0.4, "{}", m.intron_retention_index[c]);
            assert!(m.intron_retention_high[c]);
        }
        assert_eq!(m.intron_retention_high.iter().filter(|f| **f).count(), 2);
        assert_eq!(m.ir_genes_used[0], 30);
        // Global shift in every gene -> small gene dispersion.
        assert!(m.ir_gene_dispersion[3] < 0.2, "{}", m.ir_gene_dispersion[3]);
    }

    #[test]
    fn too_few_genes_is_undefined() {
        let l = layers(30, 10, 30, 10, &[]);
        let m = compute(&l, &Strata::global(30));
        assert_eq!(m.undefined_cells, 30);
        assert!(m.intron_retention_index.iter().all(|v| v.is_nan()));
        assert_eq!(m.ir_genes_used[0], 10);
    }

    #[test]
    fn reference_ratio_is_pooled_and_needs_enough_cells() {
        let l = layers(60, 25, 100, 50, &[]);
        let members: Vec<usize> = (0..60).collect();
        let r = gene_reference(&l, &members, 25);
        for g in 0..25 {
            assert!((r.p[g] - 1.0 / 3.0).abs() < 0.02, "{}", r.p[g]);
        }
        // Fewer than MIN_CELLS_PER_GENE members -> undefined reference.
        let few: Vec<usize> = (0..5).collect();
        let r = gene_reference(&l, &few, 25);
        assert!(r.p.iter().all(|p| p.is_nan()));
    }

    #[test]
    fn shrinkage_is_modest_for_well_covered_genes() {
        // Cell with 50 UMIs at ratio 0.4 against a reference of 0.25 with
        // K = 10 keeps most of its own signal: (20 + 2.5) / (50 + 10) = 0.375.
        let p = 0.25;
        let ir = (20.0 + PRIOR_STRENGTH * p) / (50.0 + PRIOR_STRENGTH);
        assert!((ir - 0.375).abs() < 1e-9);
    }
}
