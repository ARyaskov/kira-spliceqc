use crate::model::cell_cycle::{CellCycleMetrics, CellCyclePhase};
use crate::model::cell_qc::CellQc;
use crate::model::collapse::{SpliceosomeCollapseMetrics, SpliceosomeCollapseStatus};
use crate::model::cryptic_risk::CrypticSplicingRiskMetrics;
use crate::model::intron_retention::IntronRetentionMetrics;
use crate::model::junctions::JunctionMetrics;
use crate::model::sis::{SpliceIntegrityClass, SpliceIntegrityMetrics};
use crate::model::unspliced::UnsplicedMetrics;
use crate::reference::Strata;
use crate::stats::robust::median;

#[allow(clippy::too_many_arguments)]
pub fn format_summary(
    metrics: &SpliceIntegrityMetrics,
    cryptic: Option<&CrypticSplicingRiskMetrics>,
    collapse: Option<&SpliceosomeCollapseMetrics>,
    cell_cycle: &CellCycleMetrics,
    cell_qc: &CellQc,
    unspliced: Option<&UnsplicedMetrics>,
    intron_retention: Option<&IntronRetentionMetrics>,
    junctions: Option<&JunctionMetrics>,
    strata: &Strata,
    experimental: bool,
) -> String {
    let n_cells = metrics.sis.len();
    let cycling_line = {
        let known = cell_cycle
            .phase
            .iter()
            .filter(|p| **p != CellCyclePhase::Unknown)
            .count();
        if known == 0 {
            "Cycling (S/G2M, Tirosh 2016): N/A (cell-cycle genes not mapped)\n".to_string()
        } else {
            let cycling = cell_cycle.cycling.iter().filter(|c| **c).count();
            format!(
                "Cycling (S/G2M, Tirosh 2016): {} ({:.1}%)\n",
                cycling,
                100.0 * cycling as f32 / known as f32
            )
        }
    };
    let qc_line = format!(
        "Cell QC: {} low-depth (< {} UMIs or < {} genes), {} doublets; excluded from reference norms\n",
        cell_qc.n_low_depth(),
        cell_qc.min_counts,
        cell_qc.min_genes,
        cell_qc.n_doublet()
    );
    let reference = match &strata.column {
        Some(col) => format!(
            "{qc_line}Reference: {} by {} ({} strata, {} cells folded into global)\n",
            strata.mode.as_str(),
            col,
            strata.n_strata() - 1,
            strata.folded_cells
        ),
        None => format!(
            "{qc_line}Reference: {} (no stratification column)\n",
            strata.mode.as_str()
        ),
    };
    let tier_b = match junctions {
        Some(j) => format!(
            "Junctions (L2): {} junctions, {} annotated, {} cryptic-3'SS, {} skip; undefined in {} cells; cryptic-high {}, skip-high {}, shift-high {}\n",
            j.n_junctions,
            j.n_annotated,
            j.n_cryptic_acceptor_junctions,
            j.n_skip_junctions,
            j.undefined_cells,
            j.cryptic_3ss_high.iter().filter(|f| **f).count(),
            j.exon_skip_high.iter().filter(|f| **f).count(),
            j.splice_site_shift_high.iter().filter(|f| **f).count()
        ),
        None => String::new(),
    };
    let tier_a = match unspliced {
        Some(u) => {
            let med = median(&u.unspliced_fraction);
            let flagged = u.nuclear_fraction_flag.iter().filter(|f| **f).count();
            let iri = match intron_retention {
                Some(ir) => format!(
                    "Intron retention index: median {:.3}, undefined in {} cells, high flags: {}\n",
                    median(&ir.intron_retention_index),
                    ir.undefined_cells,
                    ir.intron_retention_high.iter().filter(|f| **f).count()
                ),
                None => String::new(),
            };
            format!(
                "Input levels: L0, L1 ({})\n{}Unspliced fraction: median {:.3}, undefined in {} cells, nuclear-fraction flags: {}\n{}{}",
                u.source, reference, med, u.undefined_cells, flagged, iri, tier_b
            )
        }
        None => format!(
            "Input levels: L0 (no spliced/unspliced layers; Tier A metrics unavailable)\n{}{}",
            reference, tier_b
        ),
    };
    if !experimental {
        let undefined = metrics.sis.iter().filter(|v| !v.is_finite()).count();
        return format!(
            "kira-spliceqc summary\n---------------------\nCells analyzed: {}\n{}{}Cells with undefined expression signatures: {}\n\nComposite signatures (SIS classes, SOS/RLR/SII, cryptic risk, collapse) are\nexperimental and not written; pass --experimental-signatures to include them.\n",
            n_cells, tier_a, cycling_line, undefined
        );
    }
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

    format!(
        "kira-spliceqc summary\n---------------------\nCells analyzed: {}\n{}{}\nIntegrity classes:\n  Intact:    {:>4} ({:.1}%)\n  Stressed:  {:>4} ({:.1}%)\n  Impaired:  {:>4} ({:.1}%)\n  Broken:    {:>4} ({:.1}%)\n\nMedian SIS: {:.2}\nFailure fraction (Impaired+Broken): {:.1}%\n\nCryptic splicing risk > 0.7: {}\nSpliceosome collapse: {}\n",
        n_cells,
        tier_a,
        cycling_line,
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
        fmt_pct(collapse_pct)
    )
}

fn fmt_pct(value: Option<Option<f32>>) -> String {
    match value {
        Some(Some(pct)) => format!("{pct:.1}%"),
        _ => "N/A".to_string(),
    }
}
