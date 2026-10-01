//! Tier A: per-cell unspliced fraction from spliced/unspliced layers.
//!
//! `UF_c = U_c / (S_c + U_c)` with a Wilson 95 % confidence interval. This
//! is a direct observation of immature transcripts per cell (La Manno et al.
//! 2018), not an expression signature. Its expected level depends on the
//! protocol (nuclei vs whole cells) and cell type, so interpretation is
//! relative to a reference stratum (see the reference module).

use rayon::prelude::*;

use crate::expression::SplicedUnspliced;
use crate::model::unspliced::UnsplicedMetrics;
use crate::reference::external::ReferenceFile;
use crate::reference::{
    Strata, apply_proportion_norms, flag_outliers, logit_deviation_by_stratum, proportion_norms,
};

/// Cells with fewer spliced + unspliced UMIs get an undefined fraction: a
/// binomial proportion on fewer trials has a Wilson interval wider than
/// ~0.2 and carries no usable information.
pub const MIN_LAYER_UMIS: u64 = 100;

/// z for a 95 % two-sided interval.
const Z95: f64 = 1.959_963_984_540_054;

/// With an external reference the deviations use the file's per-stratum
/// norms (cells are already assigned to the reference strata); the
/// dataset's own norms are still computed and reported.
pub fn compute(
    layers: &SplicedUnspliced,
    strata: &Strata,
    external: Option<&ReferenceFile>,
) -> UnsplicedMetrics {
    let n_cells = layers.n_cells();
    debug_assert_eq!(strata.n_cells(), n_cells);

    let rows: Vec<(u64, u64, u64, f32, f32, f32)> = (0..n_cells)
        .into_par_iter()
        .map(|cell| {
            let s = layers.spliced.cell_total(cell);
            let u = layers.unspliced.cell_total(cell);
            let a = layers.ambiguous.as_ref().map_or(0, |m| m.cell_total(cell));
            let n = s + u;
            if n < MIN_LAYER_UMIS {
                (s, u, a, f32::NAN, f32::NAN, f32::NAN)
            } else {
                let (lo, hi) = wilson_interval(u, n, Z95);
                (s, u, a, (u as f64 / n as f64) as f32, lo as f32, hi as f32)
            }
        })
        .collect();

    let mut spliced_umis = Vec::with_capacity(n_cells);
    let mut unspliced_umis = Vec::with_capacity(n_cells);
    let mut ambiguous_umis = Vec::with_capacity(n_cells);
    let mut unspliced_fraction = Vec::with_capacity(n_cells);
    let mut ci_low = Vec::with_capacity(n_cells);
    let mut ci_high = Vec::with_capacity(n_cells);
    let mut undefined_cells = 0usize;
    for (s, u, a, f, lo, hi) in rows {
        if !f.is_finite() {
            undefined_cells += 1;
        }
        spliced_umis.push(s);
        unspliced_umis.push(u);
        ambiguous_umis.push(a);
        unspliced_fraction.push(f);
        ci_low.push(lo);
        ci_high.push(hi);
    }

    let trials: Vec<u64> = spliced_umis
        .iter()
        .zip(&unspliced_umis)
        .map(|(s, u)| s + u)
        .collect();
    let (_, reference) = logit_deviation_by_stratum(&unspliced_fraction, &trials, strata);
    let norms = proportion_norms(&unspliced_fraction, &trials, strata);
    let (mut unspliced_fraction_dev, norm_source) = match external {
        Some(file) => (
            apply_proportion_norms(
                &unspliced_fraction,
                &trials,
                &strata.labels,
                &file.unspliced_norms(),
            ),
            "external",
        ),
        None => {
            let own: Vec<_> = norms.iter().map(|n| Some(*n)).collect();
            (
                apply_proportion_norms(&unspliced_fraction, &trials, &strata.labels, &own),
                "internal",
            )
        }
    };
    strata.blank_excluded(&mut unspliced_fraction_dev);
    let nuclear_fraction_flag = flag_outliers(&unspliced_fraction_dev, strata, -1.0);

    UnsplicedMetrics {
        source: layers.source.clone(),
        min_layer_umis: MIN_LAYER_UMIS,
        has_ambiguous: layers.ambiguous.is_some(),
        spliced_umis,
        unspliced_umis,
        ambiguous_umis,
        unspliced_fraction,
        unspliced_fraction_ci_low: ci_low,
        unspliced_fraction_ci_high: ci_high,
        cells_without_layers: layers.cells_without_layers,
        undefined_cells,
        unspliced_fraction_dev,
        nuclear_fraction_flag,
        reference,
        norms,
        norm_source,
    }
}

/// Wilson score interval for `successes / trials` at normal quantile `z`.
pub fn wilson_interval(successes: u64, trials: u64, z: f64) -> (f64, f64) {
    if trials == 0 {
        return (f64::NAN, f64::NAN);
    }
    let n = trials as f64;
    let p = successes as f64 / n;
    let z2 = z * z;
    let denom = 1.0 + z2 / n;
    let center = (p + z2 / (2.0 * n)) / denom;
    let half = z * (p * (1.0 - p) / n + z2 / (4.0 * n * n)).sqrt() / denom;
    // Exact bounds at the extremes (the algebra leaves ~1e-18 residues).
    let lo = if successes == 0 {
        0.0
    } else {
        (center - half).max(0.0)
    };
    let hi = if successes == trials {
        1.0
    } else {
        (center + half).min(1.0)
    };
    (lo, hi)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::expression::LayerMatrix;

    #[test]
    fn wilson_matches_reference_values() {
        // 20 of 100: Wilson 95 % = [0.1334, 0.2888].
        let (lo, hi) = wilson_interval(20, 100, Z95);
        assert!((lo - 0.1334).abs() < 1e-3, "{lo}");
        assert!((hi - 0.2888).abs() < 1e-3, "{hi}");
        // Boundaries: 0 successes -> lower bound exactly 0.
        let (lo, hi) = wilson_interval(0, 50, Z95);
        assert_eq!(lo, 0.0);
        assert!(hi > 0.0 && hi < 0.1);
        let (lo, hi) = wilson_interval(50, 50, Z95);
        assert!(lo > 0.9);
        assert_eq!(hi, 1.0);
        assert!(wilson_interval(0, 0, Z95).0.is_nan());
    }

    #[test]
    fn fraction_and_threshold() {
        let spliced = LayerMatrix::from_triplets(2, 2, vec![(0, 0, 150), (0, 1, 10)]);
        let unspliced = LayerMatrix::from_triplets(2, 2, vec![(1, 0, 50), (1, 1, 10)]);
        let ambiguous = LayerMatrix::from_triplets(2, 2, vec![(1, 0, 8)]);
        let layers = SplicedUnspliced {
            spliced,
            unspliced,
            ambiguous: Some(ambiguous),
            source: "test".to_string(),
            cells_without_layers: 0,
        };
        let m = compute(&layers, &Strata::global(2), None);
        assert_eq!(m.spliced_umis, vec![150, 10]);
        assert_eq!(m.unspliced_umis, vec![50, 10]);
        assert_eq!(m.ambiguous_umis, vec![8, 0]);
        assert!((m.unspliced_fraction[0] - 0.25).abs() < 1e-6);
        assert!(m.unspliced_fraction_ci_low[0] < 0.25 && m.unspliced_fraction_ci_high[0] > 0.25);
        // 20 UMIs < MIN_LAYER_UMIS -> undefined.
        assert!(m.unspliced_fraction[1].is_nan());
        assert_eq!(m.undefined_cells, 1);
        assert!(m.has_ambiguous);
    }
}
