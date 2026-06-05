use crate::expression::ExpressionMatrix;
use crate::stats::robust;

use super::panels::PANEL_TRIM_FRAC;

pub use crate::stats::robust::RobustRef;

/// Trimmed mean of `log_cp10k` over a panel of sorted gene indices.
///
/// Reuses the caller-supplied `scratch` to avoid per-cell allocations.
/// Uses one sparse merge pass via `ExpressionMatrix::gather_panel_log_cp10k`
/// instead of `panel_indices.len()` binary searches.
pub fn panel_trimmed_mean(
    matrix: &dyn ExpressionMatrix,
    panel_indices: &[u32],
    cell: usize,
    min_genes: usize,
    scratch: &mut Vec<f32>,
) -> f32 {
    scratch.clear();
    if panel_indices.is_empty() {
        return f32::NAN;
    }
    scratch.reserve(panel_indices.len());
    matrix.gather_panel_log_cp10k(panel_indices, cell, scratch);

    let n = scratch.len();
    if n < min_genes {
        return f32::NAN;
    }

    let trim = ((n as f32) * PANEL_TRIM_FRAC).floor() as usize;
    let start = trim.min(n);
    let end = n.saturating_sub(trim).max(start);
    if start >= end {
        return f32::NAN;
    }

    // O(n) two-sided partition: place the lower `start` smallest values into [0..start)
    // and the upper `n-end` largest into [end..n). The middle slice is then the trim window.
    if start > 0 {
        scratch.select_nth_unstable_by(start - 1, |a, b| a.total_cmp(b));
    }
    if end < n {
        // After the first partition only the right half [start..n) is unsorted relative to
        // the upper boundary, so the kth-from-start position is `end-1-start`.
        let kth = end - 1 - start;
        scratch[start..].select_nth_unstable_by(kth, |a, b| a.total_cmp(b));
    }

    let slice = &scratch[start..end];
    let sum: f32 = slice.iter().copied().sum();
    sum / slice.len() as f32
}

#[inline]
pub fn robust_z(values: &[f32]) -> (Vec<f32>, RobustRef) {
    robust::robust_z(values)
}

#[inline]
pub fn percentile(values: &[f32], q: f32) -> f32 {
    robust::percentile(values, q)
}
