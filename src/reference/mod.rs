//! Reference strata and deviation statistics.
//!
//! Every metric is compared with the reference of its own stratum, never with
//! the median of the whole dataset unless nothing better exists. Strata come
//! from cell metadata (`cell_type`, then `cluster`), fall back to a single
//! global stratum, and small strata (< `MIN_STRATUM_CELLS`) fold into the
//! global one. `summary.json` records which mode was used.

use std::collections::BTreeMap;

use tracing::{info, warn};

use crate::input::metadata::{CELL_TYPE_ALIASES, CLUSTER_ALIASES, CellMetadata};
use crate::stats::robust::{mad, median};

/// Strata smaller than this are folded into the global stratum.
pub const MIN_STRATUM_CELLS: usize = 50;
/// Name of the fold-back stratum.
pub const GLOBAL_STRATUM: &str = "global";
/// |deviation| at or above which a cell is a candidate outlier.
pub const DEVIATION_THRESHOLD: f32 = 3.0;
/// False-discovery rate for the Benjamini-Hochberg step on outlier flags.
pub const FLAG_FDR: f64 = 0.05;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReferenceMode {
    /// Norms from a user-supplied reference file (not yet available).
    External,
    /// Norms per metadata stratum.
    Stratified,
    /// One stratum: the whole dataset. Outlier flags are still computed but
    /// classes are unreliable when a large population deviates.
    Global,
}

impl ReferenceMode {
    pub fn as_str(&self) -> &'static str {
        match self {
            ReferenceMode::External => "external",
            ReferenceMode::Stratified => "stratified",
            ReferenceMode::Global => "global",
        }
    }
}

#[derive(Debug, Clone)]
pub struct Strata {
    pub mode: ReferenceMode,
    /// Metadata column the strata came from, if any.
    pub column: Option<String>,
    /// Stratum index per cell.
    pub labels: Vec<u32>,
    /// Stratum names, indexed by label.
    pub names: Vec<String>,
    /// Cells folded into the global stratum because their own was too small.
    pub folded_cells: usize,
}

impl Strata {
    pub fn global(n_cells: usize) -> Self {
        Self {
            mode: ReferenceMode::Global,
            column: None,
            labels: vec![0; n_cells],
            names: vec![GLOBAL_STRATUM.to_string()],
            folded_cells: 0,
        }
    }

    pub fn n_strata(&self) -> usize {
        self.names.len()
    }

    pub fn n_cells(&self) -> usize {
        self.labels.len()
    }

    /// Cell indices of every stratum, in label order.
    pub fn members(&self) -> Vec<Vec<usize>> {
        let mut out = vec![Vec::new(); self.names.len()];
        for (cell, &label) in self.labels.iter().enumerate() {
            out[label as usize].push(cell);
        }
        out
    }

    pub fn sizes(&self) -> Vec<usize> {
        self.members().iter().map(|m| m.len()).collect()
    }

    /// Strata from metadata: an explicit column, else the first cell-type or
    /// cluster alias present, else global.
    pub fn from_metadata(metadata: &CellMetadata, n_cells: usize, column: Option<&str>) -> Self {
        let chosen: Option<(String, Vec<String>)> = match column {
            Some(name) => match metadata.resolve(&[name]) {
                Some((found, values)) => Some((found.to_string(), values.to_vec())),
                None => {
                    warn!(column = name, "stratification column not found in metadata; using global reference");
                    None
                }
            },
            None => metadata
                .resolve(CELL_TYPE_ALIASES)
                .or_else(|| metadata.resolve(CLUSTER_ALIASES))
                .map(|(name, values)| (name.to_string(), values.to_vec())),
        };
        let Some((name, values)) = chosen else {
            info!("no stratification column; global reference");
            return Self::global(n_cells);
        };
        if values.len() != n_cells {
            warn!(column = name.as_str(), "stratification column length mismatch; using global reference");
            return Self::global(n_cells);
        }

        // Count per value; empty values and small groups fold into global.
        let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
        for v in &values {
            *counts.entry(v.as_str()).or_insert(0) += 1;
        }
        let mut names: Vec<String> = vec![GLOBAL_STRATUM.to_string()];
        let mut index: BTreeMap<&str, u32> = BTreeMap::new();
        for (value, count) in &counts {
            if !value.is_empty() && *count >= MIN_STRATUM_CELLS {
                index.insert(value, names.len() as u32);
                names.push((*value).to_string());
            }
        }
        let mut folded = 0usize;
        let labels: Vec<u32> = values
            .iter()
            .map(|v| match index.get(v.as_str()) {
                Some(&l) => l,
                None => {
                    folded += 1;
                    0
                }
            })
            .collect();
        if names.len() == 1 {
            warn!(
                column = name.as_str(),
                min_cells = MIN_STRATUM_CELLS,
                "every stratum is smaller than the minimum; using global reference"
            );
            return Self::global(n_cells);
        }
        // Drop the fold-back stratum from the name list when unused? Keep it:
        // label 0 is always "global" so downstream code stays simple.
        info!(
            column = name.as_str(),
            strata = names.len() - 1,
            folded_cells = folded,
            "stratified reference"
        );
        Self {
            mode: ReferenceMode::Stratified,
            column: Some(name),
            labels,
            names,
            folded_cells: folded,
        }
    }
}

/// Per-stratum robust reference of one metric.
#[derive(Debug, Clone)]
pub struct StratumStat {
    pub name: String,
    pub n_cells: usize,
    pub n_defined: usize,
    pub median: f32,
    pub mad: f32,
}

/// Robust z-score of `values` within each stratum: `(v - median_s) /
/// (1.4826 * MAD_s)`. NaN where the value or the stratum reference is
/// undefined, or when the stratum MAD is zero.
pub fn robust_z_by_stratum(values: &[f32], strata: &Strata) -> (Vec<f32>, Vec<StratumStat>) {
    let mut z = vec![f32::NAN; values.len()];
    let mut stats = Vec::with_capacity(strata.n_strata());
    for (label, members) in strata.members().iter().enumerate() {
        let sample: Vec<f32> = members.iter().map(|&c| values[c]).collect();
        let med = median(&sample);
        let m = mad(&sample, med);
        let n_defined = sample.iter().filter(|v| v.is_finite()).count();
        let scale = 1.4826 * m;
        if med.is_finite() && scale > 0.0 {
            for &c in members {
                if values[c].is_finite() {
                    z[c] = (values[c] - med) / scale;
                }
            }
        }
        stats.push(StratumStat {
            name: strata.names[label].clone(),
            n_cells: members.len(),
            n_defined,
            median: med,
            mad: m,
        });
    }
    (z, stats)
}

/// Deviation of a binomial proportion from its stratum on the logit scale:
/// `d = (logit p_c - logit p_s) / sqrt(1 / (n_c p_c (1 - p_c)) + tau2_s)`,
/// where `p_s` is the stratum median proportion and `tau2_s` the
/// method-of-moments overdispersion (robust variance of the cells' logits
/// minus their mean binomial variance, clamped at 0). `p_c` is clamped to
/// `[0.5 / n, 1 - 0.5 / n]` so the logit stays finite.
pub fn logit_deviation_by_stratum(
    proportions: &[f32],
    trials: &[u64],
    strata: &Strata,
) -> (Vec<f32>, Vec<StratumStat>) {
    let n_cells = proportions.len();
    let mut logits = vec![f64::NAN; n_cells];
    let mut binom_var = vec![f64::NAN; n_cells];
    for c in 0..n_cells {
        let p = proportions[c] as f64;
        let n = trials[c] as f64;
        if p.is_finite() && n > 0.0 {
            let half = 0.5 / n;
            let pc = p.clamp(half, 1.0 - half);
            logits[c] = (pc / (1.0 - pc)).ln();
            binom_var[c] = 1.0 / (n * pc * (1.0 - pc));
        }
    }

    let mut d = vec![f32::NAN; n_cells];
    let mut stats = Vec::with_capacity(strata.n_strata());
    for (label, members) in strata.members().iter().enumerate() {
        let props: Vec<f32> = members.iter().map(|&c| proportions[c]).collect();
        let med_p = median(&props);
        let mad_p = mad(&props, med_p);
        let member_logits: Vec<f32> = members.iter().map(|&c| logits[c] as f32).collect();
        let med_logit = median(&member_logits) as f64;
        let mad_logit = mad(&member_logits, med_logit as f32) as f64;
        let robust_var = (1.4826 * mad_logit).powi(2);
        let defined: Vec<usize> = members.iter().copied().filter(|&c| logits[c].is_finite()).collect();
        let mean_binom = if defined.is_empty() {
            f64::NAN
        } else {
            defined.iter().map(|&c| binom_var[c]).sum::<f64>() / defined.len() as f64
        };
        let tau2 = (robust_var - mean_binom).max(0.0);
        if med_logit.is_finite() && mean_binom.is_finite() {
            for &c in &defined {
                let denom = (binom_var[c] + tau2).sqrt();
                if denom > 0.0 {
                    d[c] = ((logits[c] - med_logit) / denom) as f32;
                }
            }
        }
        stats.push(StratumStat {
            name: strata.names[label].clone(),
            n_cells: members.len(),
            n_defined: defined.len(),
            median: med_p,
            mad: mad_p,
        });
    }
    (d, stats)
}

/// Two-sided normal tail probability `P(|Z| >= |d|)`.
pub fn two_sided_p(d: f64) -> f64 {
    if !d.is_finite() {
        return f64::NAN;
    }
    erfc(d.abs() / std::f64::consts::SQRT_2)
}

/// Benjamini-Hochberg adjusted p-values (NaN inputs are left NaN and do not
/// count toward `m`).
pub fn benjamini_hochberg(p: &[f64]) -> Vec<f64> {
    let mut order: Vec<usize> = (0..p.len()).filter(|&i| p[i].is_finite()).collect();
    let m = order.len();
    let mut adjusted = vec![f64::NAN; p.len()];
    if m == 0 {
        return adjusted;
    }
    order.sort_by(|&a, &b| p[a].partial_cmp(&p[b]).unwrap());
    let mut running = 1.0f64;
    for (rank_from_end, &i) in order.iter().enumerate().rev() {
        let rank = rank_from_end + 1;
        let value = (p[i] * m as f64 / rank as f64).min(1.0);
        running = running.min(value);
        adjusted[i] = running;
    }
    adjusted
}

/// Outlier flags in one direction: `sign * d >= DEVIATION_THRESHOLD` and
/// BH-adjusted two-sided p < FLAG_FDR within each stratum.
pub fn flag_outliers(deviation: &[f32], strata: &Strata, sign: f32) -> Vec<bool> {
    let mut flags = vec![false; deviation.len()];
    for members in strata.members() {
        let p: Vec<f64> = members
            .iter()
            .map(|&c| two_sided_p(deviation[c] as f64))
            .collect();
        let adj = benjamini_hochberg(&p);
        for (k, &c) in members.iter().enumerate() {
            let d = deviation[c];
            flags[c] = d.is_finite() && sign * d >= DEVIATION_THRESHOLD && adj[k] < FLAG_FDR;
        }
    }
    flags
}

/// Complementary error function (Abramowitz-Stegun 7.1.26 refined with one
/// continued-fraction step; absolute error < 1e-9 over the working range).
fn erfc(x: f64) -> f64 {
    // Numerical Recipes erfcc: fractional error < 1.2e-7 everywhere.
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let r = t * (-z * z - 1.265_512_23
        + t * (1.000_023_68
            + t * (0.374_091_96
                + t * (0.096_784_18
                    + t * (-0.186_288_06
                        + t * (0.278_868_07
                            + t * (-1.135_203_98
                                + t * (1.488_515_87 + t * (-0.822_152_23 + t * 0.170_872_77)))))))))
        .exp();
    if x >= 0.0 { r } else { 2.0 - r }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn metadata(col: &str, values: &[&str]) -> CellMetadata {
        let mut m = CellMetadata::default();
        m.columns
            .insert(col.to_string(), values.iter().map(|s| s.to_string()).collect());
        m
    }

    #[test]
    fn small_strata_fold_into_global() {
        let mut values: Vec<&str> = vec!["T"; 60];
        values.extend(vec!["B"; 10]);
        values.extend(vec![""; 5]);
        let s = Strata::from_metadata(&metadata("cell_type", &values), 75, None);
        assert_eq!(s.mode, ReferenceMode::Stratified);
        assert_eq!(s.column.as_deref(), Some("cell_type"));
        assert_eq!(s.names, vec!["global", "T"]);
        assert_eq!(s.folded_cells, 15);
        assert_eq!(s.sizes(), vec![15, 60]);
    }

    #[test]
    fn global_when_no_column_or_all_small() {
        let s = Strata::from_metadata(&CellMetadata::default(), 10, None);
        assert_eq!(s.mode, ReferenceMode::Global);
        let s = Strata::from_metadata(&metadata("cluster", &["a"; 10]), 10, None);
        assert_eq!(s.mode, ReferenceMode::Global);
        let s = Strata::from_metadata(&metadata("cluster", &["a"; 60]), 60, Some("missing"));
        assert_eq!(s.mode, ReferenceMode::Global);
    }

    #[test]
    fn robust_z_is_per_stratum() {
        // Stratum A around 1.0, stratum B around 10.0, each with spread
        // {-0.5, 0, +0.5} so the MAD is 0.5 and the median cell is exact.
        let values: Vec<f32> = (0..120)
            .map(|i| {
                let base = if i < 60 { 1.0 } else { 10.0 };
                base + ((i % 3) as f32 - 1.0) * 0.5
            })
            .collect();
        let mut labels: Vec<&str> = vec!["A"; 60];
        labels.extend(vec!["B"; 60]);
        let s = Strata::from_metadata(&metadata("cell_type", &labels), 120, None);
        let (z, stats) = robust_z_by_stratum(&values, &s);
        assert_eq!(stats.len(), 3);
        // i % 3 == 1 -> at the stratum median -> z = 0 in both strata.
        assert_eq!(z[1], 0.0);
        assert_eq!(z[61], 0.0);
        // +0.5 above the median with MAD 0.5 -> z = 1 / 1.4826.
        assert!((z[2] - 1.0 / 1.4826).abs() < 1e-4);
        assert!((z[62] - 1.0 / 1.4826).abs() < 1e-4);
        // A stratum-A value would be a huge outlier under a global reference
        // but is ordinary within its stratum.
        assert!(z.iter().all(|v| v.abs() < 2.0));
    }

    #[test]
    fn logit_deviation_flags_low_outlier_only() {
        // 100 cells at p = 0.3 (n = 400), one damaged cell at p = 0.02.
        let mut p = vec![0.3f32; 101];
        let n = vec![400u64; 101];
        // add spread so the robust variance is not degenerate
        for (i, v) in p.iter_mut().enumerate().take(100) {
            *v += ((i % 7) as f32 - 3.0) * 0.01;
        }
        p[100] = 0.02;
        let s = Strata::global(101);
        let (d, stats) = logit_deviation_by_stratum(&p, &n, &s);
        assert_eq!(stats[0].n_defined, 101);
        assert!(d[100] < -3.0, "{}", d[100]);
        let low = flag_outliers(&d, &s, -1.0);
        assert!(low[100]);
        assert_eq!(low.iter().filter(|f| **f).count(), 1);
        let high = flag_outliers(&d, &s, 1.0);
        assert!(!high.iter().any(|f| *f));
    }

    #[test]
    fn bh_and_normal_tail() {
        assert!((two_sided_p(1.959_964) - 0.05).abs() < 1e-4);
        assert!((two_sided_p(0.0) - 1.0).abs() < 1e-6);
        // m = 4 finite p-values: ranks 0.01 (1), 0.03 (2), 0.04 (3), 0.5 (4).
        // adj = min over higher ranks of p * m / rank:
        // 0.5 -> 0.5; 0.04 -> 0.0533; 0.03 -> min(0.06, 0.0533); 0.01 -> 0.04.
        let adj = benjamini_hochberg(&[0.01, 0.04, 0.03, f64::NAN, 0.5]);
        assert!((adj[0] - 0.04).abs() < 1e-9);
        assert!((adj[1] - 0.04 * 4.0 / 3.0).abs() < 1e-9);
        assert!((adj[2] - 0.04 * 4.0 / 3.0).abs() < 1e-9);
        assert!(adj[3].is_nan());
        assert!((adj[4] - 0.5).abs() < 1e-9);
    }
}
