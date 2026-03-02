use std::cmp::Ordering;

use crate::expression::ExpressionMatrix;

use super::panels::{EPS_ROBUST, PANEL_TRIM_FRAC};

#[derive(Debug, Clone, Copy)]
pub struct RobustRef {
    pub median: f32,
    pub mad: f32,
}

pub fn panel_trimmed_mean(
    matrix: &dyn ExpressionMatrix,
    panel_indices: &[u32],
    cell: usize,
    min_genes: usize,
    scratch: &mut Vec<f32>,
) -> f32 {
    scratch.clear();
    scratch.reserve(panel_indices.len());
    for &gene_idx in panel_indices {
        scratch.push(matrix.log_cp10k(gene_idx as usize, cell));
    }

    if scratch.len() < min_genes {
        return f32::NAN;
    }

    scratch.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let n = scratch.len();
    let trim = ((n as f32) * PANEL_TRIM_FRAC).floor() as usize;
    let start = trim.min(n);
    let end = n.saturating_sub(trim).max(start);
    let slice = &scratch[start..end];
    if slice.is_empty() {
        return f32::NAN;
    }
    let sum: f32 = slice.iter().copied().sum();
    sum / slice.len() as f32
}

pub fn robust_z(values: &[f32]) -> (Vec<f32>, RobustRef) {
    let finite = values
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .collect::<Vec<_>>();
    if finite.is_empty() {
        return (
            vec![f32::NAN; values.len()],
            RobustRef {
                median: f32::NAN,
                mad: f32::NAN,
            },
        );
    }

    let median = median_sorted(finite.clone());
    let abs_dev = finite
        .iter()
        .copied()
        .map(|v| (v - median).abs())
        .collect::<Vec<_>>();
    let mad = median_sorted(abs_dev);

    if !mad.is_finite() || mad <= 0.0 {
        return (
            values
                .iter()
                .map(|v| if v.is_finite() { 0.0 } else { f32::NAN })
                .collect(),
            RobustRef { median, mad },
        );
    }

    let denom = 1.4826 * mad + EPS_ROBUST;
    let z = values
        .iter()
        .map(|v| {
            if v.is_finite() {
                (*v - median) / denom
            } else {
                f32::NAN
            }
        })
        .collect();

    (z, RobustRef { median, mad })
}

pub fn percentile(values: &[f32], q: f32) -> f32 {
    let mut finite = values
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .collect::<Vec<_>>();
    if finite.is_empty() {
        return f32::NAN;
    }
    finite.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let max_idx = finite.len() - 1;
    let pos = ((max_idx as f32) * q).round() as usize;
    finite[pos.min(max_idx)]
}

fn median_sorted(mut values: Vec<f32>) -> f32 {
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let mid = values.len() / 2;
    if values.len() % 2 == 1 {
        values[mid]
    } else {
        0.5 * (values[mid - 1] + values[mid])
    }
}
