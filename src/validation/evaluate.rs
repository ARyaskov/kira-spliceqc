//! Scores a run's per-cell outputs against a truth table.
//!
//! For every (truth column, metric, flag) pair: AUROC (Mann-Whitney with
//! average ranks), AUPRC (precision at every positive), and for the flag
//! precision, recall, F1 and the false-positive rate among negatives.
//! Cells with an undefined metric are excluded from AUROC/AUPRC and counted.
//! With a `cell_type` column in the truth table the same numbers are also
//! reported per stratum.

use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::path::Path;

use serde::Serialize;

use crate::input::error::InputError;

/// One scoring pair: `truth` is a boolean column of the truth table,
/// `metric` a numeric column of cells.tsv (higher = more positive after
/// `sign`), `flag` an optional boolean column of cells.tsv.
#[derive(Debug, Clone, Serialize)]
pub struct Pair {
    pub truth: String,
    pub metric: String,
    pub flag: Option<String>,
    pub sign: f64,
}

impl Pair {
    /// `truth:metric[:flag[:sign]]`.
    pub fn parse(spec: &str) -> Result<Self, InputError> {
        let parts: Vec<&str> = spec.split(':').collect();
        if parts.len() < 2 || parts.len() > 4 {
            return Err(InputError::UnsupportedInput(format!("pair spec {spec:?}: expected truth:metric[:flag[:sign]]")));
        }
        let sign = match parts.get(3) {
            Some(s) => s.parse::<f64>().map_err(|_| InputError::UnsupportedInput(format!("pair spec {spec:?}: bad sign")))?,
            None => 1.0,
        };
        Ok(Self {
            truth: parts[0].to_string(),
            metric: parts[1].to_string(),
            flag: parts.get(2).filter(|f| !f.is_empty()).map(|f| f.to_string()),
            sign,
        })
    }
}

/// Default pairs for the `simulate` truth columns.
pub fn default_pairs() -> Vec<Pair> {
    [
        ("truth_cryptic", "cryptic_3ss_fraction_dev", "cryptic_3ss_high", 1.0),
        ("truth_ir", "intron_retention_index_dev", "intron_retention_high", 1.0),
        ("truth_damaged", "unspliced_fraction_dev", "nuclear_fraction_flag", -1.0),
        ("truth_skip", "exon_skip_fraction_dev", "exon_skip_high", 1.0),
    ]
    .iter()
    .map(|(t, m, f, s)| Pair {
        truth: t.to_string(),
        metric: m.to_string(),
        flag: Some(f.to_string()),
        sign: *s,
    })
    .collect()
}

#[derive(Debug, Clone, Serialize)]
pub struct PairResult {
    pub truth: String,
    pub metric: String,
    pub flag: Option<String>,
    pub stratum: String,
    pub n_positive: usize,
    pub n_negative: usize,
    pub n_undefined: usize,
    pub auroc: Option<f64>,
    pub auprc: Option<f64>,
    pub flag_precision: Option<f64>,
    pub flag_recall: Option<f64>,
    pub flag_f1: Option<f64>,
    pub flag_fpr: Option<f64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ValidationReport {
    pub tool_version: &'static str,
    pub cells_file: String,
    pub truth_file: String,
    pub n_cells_scored: usize,
    pub results: Vec<PairResult>,
}

struct Table {
    header: Vec<String>,
    rows: Vec<Vec<String>>,
}

fn read_table(path: &Path) -> Result<Table, InputError> {
    let text = fs::read_to_string(path).map_err(|e| InputError::io(path, e))?;
    let mut lines = text.lines();
    let header: Vec<String> = lines
        .next()
        .ok_or_else(|| InputError::UnsupportedInput(format!("{}: empty", path.display())))?
        .split('\t')
        .map(str::to_string)
        .collect();
    let rows = lines
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.split('\t').map(str::to_string).collect())
        .collect();
    Ok(Table { header, rows })
}

fn truthy(v: &str) -> bool {
    matches!(v.trim().to_ascii_lowercase().as_str(), "true" | "1" | "yes")
}

fn auroc(scores: &[f64], labels: &[bool]) -> Option<f64> {
    let n_pos = labels.iter().filter(|l| **l).count();
    let n_neg = labels.len() - n_pos;
    if n_pos == 0 || n_neg == 0 {
        return None;
    }
    let mut order: Vec<usize> = (0..scores.len()).collect();
    order.sort_by(|&a, &b| scores[a].partial_cmp(&scores[b]).unwrap());
    let mut rank_sum_pos = 0.0;
    let mut i = 0;
    while i < order.len() {
        let mut j = i;
        while j + 1 < order.len() && scores[order[j + 1]] == scores[order[i]] {
            j += 1;
        }
        let avg = (i + j) as f64 / 2.0 + 1.0;
        for &idx in &order[i..=j] {
            if labels[idx] {
                rank_sum_pos += avg;
            }
        }
        i = j + 1;
    }
    Some((rank_sum_pos - (n_pos * (n_pos + 1)) as f64 / 2.0) / (n_pos * n_neg) as f64)
}

fn auprc(scores: &[f64], labels: &[bool]) -> Option<f64> {
    let n_pos = labels.iter().filter(|l| **l).count();
    if n_pos == 0 || n_pos == labels.len() {
        return None;
    }
    let mut order: Vec<usize> = (0..scores.len()).collect();
    order.sort_by(|&a, &b| scores[b].partial_cmp(&scores[a]).unwrap());
    let mut tp = 0usize;
    let mut sum_precision = 0.0;
    for (k, &idx) in order.iter().enumerate() {
        if labels[idx] {
            tp += 1;
            sum_precision += tp as f64 / (k + 1) as f64;
        }
    }
    Some(sum_precision / n_pos as f64)
}

fn score_pair(pair: &Pair, stratum: &str, rows: &[(f64, bool, Option<bool>)], n_undefined: usize) -> PairResult {
    let scores: Vec<f64> = rows.iter().map(|r| r.0 * pair.sign).collect();
    let labels: Vec<bool> = rows.iter().map(|r| r.1).collect();
    let n_positive = labels.iter().filter(|l| **l).count();
    let n_negative = labels.len() - n_positive;
    let (mut tp, mut fp, mut fnn, mut tn) = (0usize, 0usize, 0usize, 0usize);
    let mut has_flag = false;
    for (_, truth, flag) in rows {
        if let Some(f) = flag {
            has_flag = true;
            match (f, truth) {
                (true, true) => tp += 1,
                (true, false) => fp += 1,
                (false, true) => fnn += 1,
                (false, false) => tn += 1,
            }
        }
    }
    let ratio = |a: usize, b: usize| if b > 0 { Some(a as f64 / b as f64) } else { None };
    let precision = if has_flag { ratio(tp, tp + fp) } else { None };
    let recall = if has_flag { ratio(tp, tp + fnn) } else { None };
    let f1 = match (precision, recall) {
        (Some(p), Some(r)) if p + r > 0.0 => Some(2.0 * p * r / (p + r)),
        _ => None,
    };
    PairResult {
        truth: pair.truth.clone(),
        metric: pair.metric.clone(),
        flag: pair.flag.clone(),
        stratum: stratum.to_string(),
        n_positive,
        n_negative,
        n_undefined,
        auroc: auroc(&scores, &labels),
        auprc: auprc(&scores, &labels),
        flag_precision: precision,
        flag_recall: recall,
        flag_f1: f1,
        flag_fpr: if has_flag { ratio(fp, fp + tn) } else { None },
    }
}

pub fn evaluate(cells_tsv: &Path, truth_tsv: &Path, pairs: &[Pair]) -> Result<ValidationReport, InputError> {
    let cells = read_table(cells_tsv)?;
    let truth = read_table(truth_tsv)?;
    let col = |t: &Table, name: &str| t.header.iter().position(|h| h == name);
    let cell_name = col(&cells, "cell_name")
        .ok_or_else(|| InputError::UnsupportedInput("cells.tsv has no cell_name column".to_string()))?;
    let truth_barcode = col(&truth, "barcode")
        .ok_or_else(|| InputError::UnsupportedInput("truth table has no barcode column".to_string()))?;
    let truth_stratum = col(&truth, "cell_type");
    let truth_index: HashMap<&str, &Vec<String>> = truth.rows.iter().map(|r| (r[truth_barcode].as_str(), r)).collect();

    let mut results = Vec::new();
    let mut n_scored = 0usize;
    for pair in pairs {
        let (Some(ti), Some(mi)) = (col(&truth, &pair.truth), col(&cells, &pair.metric)) else {
            tracing::warn!(truth = pair.truth.as_str(), metric = pair.metric.as_str(), "pair skipped: column missing");
            continue;
        };
        let fi = pair.flag.as_deref().and_then(|f| col(&cells, f));
        // stratum -> rows; "all" collects everything
        let mut by_stratum: BTreeMap<String, Vec<(f64, bool, Option<bool>)>> = BTreeMap::new();
        let mut undefined: BTreeMap<String, usize> = BTreeMap::new();
        for row in &cells.rows {
            let Some(t) = truth_index.get(row[cell_name].as_str()) else { continue };
            n_scored += 1;
            let label = truthy(&t[ti]);
            let flag = fi.map(|i| truthy(&row[i]));
            let stratum = truth_stratum.map(|i| t[i].clone()).unwrap_or_default();
            let keys: Vec<String> = if stratum.is_empty() { vec!["all".to_string()] } else { vec!["all".to_string(), stratum] };
            match row[mi].parse::<f64>() {
                Ok(v) if v.is_finite() => {
                    for k in keys {
                        by_stratum.entry(k).or_default().push((v, label, flag));
                    }
                }
                _ => {
                    for k in keys {
                        *undefined.entry(k).or_default() += 1;
                    }
                }
            }
        }
        for (stratum, rows) in &by_stratum {
            results.push(score_pair(pair, stratum, rows, *undefined.get(stratum).unwrap_or(&0)));
        }
    }
    Ok(ValidationReport {
        tool_version: env!("CARGO_PKG_VERSION"),
        cells_file: cells_tsv.display().to_string(),
        truth_file: truth_tsv.display().to_string(),
        n_cells_scored: n_scored / pairs.len().max(1),
        results,
    })
}

impl ValidationReport {
    pub fn to_markdown(&self) -> String {
        let fmt = |v: Option<f64>| v.map(|x| format!("{x:.3}")).unwrap_or_else(|| "-".to_string());
        let mut s = String::new();
        s.push_str(&format!(
            "# kira-spliceqc validation report\n\nTool {} · cells `{}` · truth `{}` · {} cells scored\n\n",
            self.tool_version, self.cells_file, self.truth_file, self.n_cells_scored
        ));
        s.push_str("| truth | metric | stratum | pos | neg | undefined | AUROC | AUPRC | flag | precision | recall | F1 | FPR |\n");
        s.push_str("| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: |\n");
        for r in &self.results {
            s.push_str(&format!(
                "| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |\n",
                r.truth,
                r.metric,
                r.stratum,
                r.n_positive,
                r.n_negative,
                r.n_undefined,
                fmt(r.auroc),
                fmt(r.auprc),
                r.flag.as_deref().unwrap_or("-"),
                fmt(r.flag_precision),
                fmt(r.flag_recall),
                fmt(r.flag_f1),
                fmt(r.flag_fpr),
            ));
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn auroc_and_auprc_reference_values() {
        let scores = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4];
        let labels = [true, true, false, true, false, false];
        // Positives ranked 1,2,4 of 6 -> AUC = (6+5+3 - 6)/9 = 8/9
        assert!((auroc(&scores, &labels).unwrap() - 8.0 / 9.0).abs() < 1e-9);
        // precisions at positives: 1/1, 2/2, 3/4 -> mean 0.9167
        assert!((auprc(&scores, &labels).unwrap() - (1.0 + 1.0 + 0.75) / 3.0).abs() < 1e-9);
        assert!(auroc(&scores, &[true; 6]).is_none());
        // ties get average ranks
        assert!((auroc(&[1.0, 1.0, 0.0], &[true, false, false]).unwrap() - 0.75).abs() < 1e-9);
    }

    #[test]
    fn pair_spec_parsing() {
        let p = Pair::parse("truth_x:metric_y:flag_z:-1").unwrap();
        assert_eq!(p.sign, -1.0);
        assert_eq!(p.flag.as_deref(), Some("flag_z"));
        let p = Pair::parse("a:b").unwrap();
        assert!(p.flag.is_none());
        assert!(Pair::parse("a").is_err());
        assert!(Pair::parse("a:b:c:x").is_err());
    }
}
