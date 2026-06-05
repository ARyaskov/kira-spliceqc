//! Robust statistics using O(n) quickselect (`select_nth_unstable_by`) instead
//! of O(n log n) sort. All routines treat non-finite inputs as missing.

const ROBUST_EPS: f32 = 1e-6;
const MAD_TO_SIGMA: f32 = 1.4826;

#[derive(Debug, Clone, Copy)]
pub struct RobustRef {
    pub median: f32,
    pub mad: f32,
}

/// Median of finite values. Returns NaN if no finite input.
pub fn median(values: &[f32]) -> f32 {
    let mut buf: Vec<f32> = values.iter().copied().filter(|v| v.is_finite()).collect();
    median_in_place(&mut buf)
}

/// Median Absolute Deviation around `med`. Returns NaN on empty / all-NaN.
pub fn mad(values: &[f32], med: f32) -> f32 {
    let mut buf: Vec<f32> = values
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .map(|v| (v - med).abs())
        .collect();
    median_in_place(&mut buf)
}

/// Robust z-score: `(v - median) / (1.4826 * mad + eps)`.
/// Returns the z-score vector aligned with input (NaN preserved) and the
/// `(median, mad)` reference used.
///
/// When `mad <= 0`, finite inputs map to 0.0 (collapse-safe behavior shared
/// across the stage 4/5/8-13 hot paths).
pub fn robust_z(values: &[f32]) -> (Vec<f32>, RobustRef) {
    let mut finite: Vec<f32> = values.iter().copied().filter(|v| v.is_finite()).collect();
    if finite.is_empty() {
        return (
            vec![f32::NAN; values.len()],
            RobustRef {
                median: f32::NAN,
                mad: f32::NAN,
            },
        );
    }
    let med = median_in_place(&mut finite);

    // Reuse `finite` as scratch for |v - med|.
    for v in &mut finite {
        *v = (*v - med).abs();
    }
    let m = median_in_place(&mut finite);

    if !m.is_finite() || m <= 0.0 {
        let z = values
            .iter()
            .map(|v| if v.is_finite() { 0.0 } else { f32::NAN })
            .collect();
        return (z, RobustRef { median: med, mad: m });
    }

    let denom = MAD_TO_SIGMA * m + ROBUST_EPS;
    let z = values
        .iter()
        .map(|v| {
            if v.is_finite() {
                (*v - med) / denom
            } else {
                f32::NAN
            }
        })
        .collect();
    (z, RobustRef { median: med, mad: m })
}

/// Quantile by the `round((n-1) * q)` rule on finite values.
pub fn percentile(values: &[f32], q: f32) -> f32 {
    let mut buf: Vec<f32> = values.iter().copied().filter(|v| v.is_finite()).collect();
    if buf.is_empty() {
        return f32::NAN;
    }
    let max_idx = buf.len() - 1;
    let pos = ((max_idx as f32) * q).round() as usize;
    let k = pos.min(max_idx);
    *buf.select_nth_unstable_by(k, f32_total_cmp).1
}

/// Quantile on an already-finite f64 slice.
pub fn quantile_f64(values: &mut [f64], q: f64) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let max_idx = values.len() - 1;
    let pos = ((max_idx as f64) * q).round() as usize;
    let k = pos.min(max_idx);
    *values
        .select_nth_unstable_by(k, |a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .1
}

/// In-place median; `buf` is permuted.
fn median_in_place(buf: &mut [f32]) -> f32 {
    if buf.is_empty() {
        return f32::NAN;
    }
    let n = buf.len();
    let mid = n / 2;
    let upper = *buf.select_nth_unstable_by(mid, f32_total_cmp).1;
    if n % 2 == 1 {
        upper
    } else {
        // Lower-half max via second selection (mid-1) constrained to [0, mid).
        let lower = *buf[..mid]
            .select_nth_unstable_by(mid - 1, f32_total_cmp)
            .1;
        (lower + upper) * 0.5
    }
}

#[inline(always)]
fn f32_total_cmp(a: &f32, b: &f32) -> std::cmp::Ordering {
    a.total_cmp(b)
}

/// Zero-copy view of one geneset's column in a `n_genesets × n_cells` row-major
/// activity matrix.
#[inline]
pub fn extract_geneset_slice<'a>(values: &'a [f32], geneset_idx: usize, n_cells: usize) -> &'a [f32] {
    let start = geneset_idx * n_cells;
    &values[start..start + n_cells]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn median_odd_even() {
        assert_eq!(median(&[3.0, 1.0, 2.0]), 2.0);
        assert_eq!(median(&[1.0, 2.0, 3.0, 4.0]), 2.5);
    }

    #[test]
    fn median_skips_nan() {
        assert_eq!(median(&[1.0, f32::NAN, 3.0]), 2.0);
    }

    #[test]
    fn median_empty_is_nan() {
        assert!(median(&[]).is_nan());
        assert!(median(&[f32::NAN, f32::INFINITY]).is_nan());
    }

    #[test]
    fn percentile_matches_round_rule() {
        // n=5, q=0.5 → round(4*0.5)=2 → values_sorted[2]=3.0
        assert_eq!(percentile(&[5.0, 1.0, 3.0, 2.0, 4.0], 0.5), 3.0);
        // q=0.0 → first, q=1.0 → last
        assert_eq!(percentile(&[5.0, 1.0, 3.0, 2.0, 4.0], 0.0), 1.0);
        assert_eq!(percentile(&[5.0, 1.0, 3.0, 2.0, 4.0], 1.0), 5.0);
    }

    #[test]
    fn robust_z_collapses_when_mad_zero() {
        let (z, r) = robust_z(&[1.0, 1.0, 1.0, 1.0]);
        assert_eq!(r.median, 1.0);
        assert_eq!(r.mad, 0.0);
        assert!(z.iter().all(|v| *v == 0.0));
    }

    #[test]
    fn robust_z_preserves_nan() {
        let (z, _) = robust_z(&[1.0, 2.0, f32::NAN, 3.0]);
        assert!(z[2].is_nan());
        assert!(z[0].is_finite() && z[1].is_finite() && z[3].is_finite());
    }
}
