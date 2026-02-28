use crate::model::collapse::{SpliceosomeCollapseMetrics, SpliceosomeCollapseStatus};
use crate::model::cryptic_risk::CrypticSplicingRiskMetrics;
use crate::model::sis::{SpliceIntegrityClass, SpliceIntegrityMetrics};
use crate::stats::robust::median;

pub fn format_summary(
    metrics: &SpliceIntegrityMetrics,
    cryptic: Option<&CrypticSplicingRiskMetrics>,
    collapse: Option<&SpliceosomeCollapseMetrics>,
    cell_cycle_confounded: Option<&[bool]>,
) -> String {
    let n_cells = metrics.sis.len();
    let mut intact = 0usize;
    let mut stressed = 0usize;
    let mut impaired = 0usize;
    let mut broken = 0usize;

    for class in &metrics.class {
        match class {
            SpliceIntegrityClass::Intact => intact += 1,
            SpliceIntegrityClass::Stressed => stressed += 1,
            SpliceIntegrityClass::Impaired => impaired += 1,
            SpliceIntegrityClass::Broken => broken += 1,
        }
    }

    let med = median(&metrics.sis);
    let fail = impaired + broken;
    let pct = |count: usize| {
        if n_cells == 0 {
            0.0
        } else {
            (count as f32 / n_cells as f32) * 100.0
        }
    };

    let cryptic_pct = cryptic.map(|metrics| {
        let mut flagged = 0usize;
        let mut total = 0usize;
        for value in &metrics.cryptic_risk {
            if value.is_finite() {
                total += 1;
                if *value > 0.7 {
                    flagged += 1;
                }
            }
        }
        if total == 0 {
            None
        } else {
            Some((flagged as f32 / total as f32) * 100.0)
        }
    });

    let collapse_pct = collapse.map(|metrics| {
        let mut flagged = 0usize;
        let mut total = 0usize;
        for status in &metrics.collapse_status {
            if *status != SpliceosomeCollapseStatus::Inconclusive {
                total += 1;
                if *status == SpliceosomeCollapseStatus::Collapse {
                    flagged += 1;
                }
            }
        }
        if total == 0 {
            None
        } else {
            Some((flagged as f32 / total as f32) * 100.0)
        }
    });

    let cell_cycle_pct = cell_cycle_confounded.map(|flags| {
        let mut flagged = 0usize;
        for value in flags {
            if *value {
                flagged += 1;
            }
        }
        if flags.is_empty() {
            None
        } else {
            Some((flagged as f32 / flags.len() as f32) * 100.0)
        }
    });

    format!(
        "kira-spliceqc summary\n---------------------\nCells analyzed: {}\n\nIntegrity classes:\n  Intact:    {:>4} ({:.1}%)\n  Stressed:  {:>4} ({:.1}%)\n  Impaired:  {:>4} ({:.1}%)\n  Broken:    {:>4} ({:.1}%)\n\nMedian SIS: {:.2}\nFailure fraction (Impaired+Broken): {:.1}%\n\nCryptic splicing risk > 0.7: {}\nSpliceosome collapse: {}\nCell-cycle confounded: {}\n",
        n_cells,
        intact,
        pct(intact),
        stressed,
        pct(stressed),
        impaired,
        pct(impaired),
        broken,
        pct(broken),
        med,
        pct(fail),
        fmt_pct(cryptic_pct),
        fmt_pct(collapse_pct),
        fmt_pct(cell_cycle_pct)
    )
}

fn fmt_pct(value: Option<Option<f32>>) -> String {
    match value {
        Some(Some(pct)) => format!("{pct:.1}%"),
        _ => "N/A".to_string(),
    }
}
