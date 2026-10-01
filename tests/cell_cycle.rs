//! Cell-cycle scoring: Tirosh 2016 gene lists, control-gene correction and
//! the Seurat phase rule on a synthetic matrix with known S / G2M cells.

use kira_spliceqc::expression::ExpressionMatrix;
use kira_spliceqc::genesets::controls::ControlPool;
use kira_spliceqc::metrics::cell_cycle::{G2M_GENES, S_GENES, cell_cycle_gene_ids, compute};
use kira_spliceqc::metrics::splicing_instability::panels::CONFLICT_RISK_PANEL;
use kira_spliceqc::model::cell_cycle::CellCyclePhase;

struct Dense {
    genes: Vec<String>,
    cells: Vec<String>,
    counts: Vec<Vec<u32>>, // [gene][cell]
}

impl ExpressionMatrix for Dense {
    fn n_genes(&self) -> usize {
        self.genes.len()
    }
    fn n_cells(&self) -> usize {
        self.cells.len()
    }
    fn gene_symbol(&self, g: usize) -> &str {
        &self.genes[g]
    }
    fn cell_name(&self, c: usize) -> &str {
        &self.cells[c]
    }
    fn libsize(&self, c: usize) -> u64 {
        self.counts.iter().map(|row| row[c] as u64).sum()
    }
    fn count(&self, g: usize, c: usize) -> u32 {
        self.counts[g][c]
    }
}

/// 43 S + 54 G2M + 200 filler genes; cells 0..20 are S-phase (S genes x8),
/// 20..40 are G2/M (G2M genes x8), the rest quiescent. Every gene has a
/// baseline of 2-6 counts with a deterministic jitter.
fn dataset(n_cells: usize) -> Dense {
    let mut genes: Vec<String> = S_GENES.iter().map(|s| s.to_string()).collect();
    genes.extend(G2M_GENES.iter().map(|s| s.to_string()));
    genes.extend((0..200).map(|i| format!("FILLER{i:03}")));
    let n_s = S_GENES.len();
    let n_g2m = G2M_GENES.len();
    let counts = (0..genes.len())
        .map(|g| {
            (0..n_cells)
                .map(|c| {
                    let base = 2 + ((g * 13 + c * 7) % 5) as u32;
                    let s_cell = c < 20 && g < n_s;
                    let g2m_cell = (20..40).contains(&c) && (n_s..n_s + n_g2m).contains(&g);
                    let boost = if s_cell || g2m_cell { 8 } else { 1 };
                    base * boost
                })
                .collect()
        })
        .collect();
    Dense {
        genes,
        cells: (0..n_cells).map(|c| format!("c{c}")).collect(),
        counts,
    }
}

#[test]
fn phases_follow_the_seurat_rule() {
    let m = dataset(120);
    let exclude = cell_cycle_gene_ids(&m);
    assert_eq!(exclude.len(), S_GENES.len() + G2M_GENES.len());
    let pool = ControlPool::new(m.gene_mean_log_cp10k(), &exclude);
    let cc = compute(&m, Some(&pool));
    assert_eq!(cc.s_genes_mapped, S_GENES.len());
    assert_eq!(cc.g2m_genes_mapped, G2M_GENES.len());
    for c in 0..20 {
        assert_eq!(
            cc.phase[c],
            CellCyclePhase::S,
            "cell {c}: s={} g2m={}",
            cc.s_score[c],
            cc.g2m_score[c]
        );
        assert!(cc.cycling[c]);
        assert!(cc.s_score[c] > 0.5);
    }
    for c in 20..40 {
        assert_eq!(cc.phase[c], CellCyclePhase::G2M, "cell {c}");
        assert!(cc.cycling[c]);
    }
    // Quiescent cells: both scores near the control background. (Boosting
    // the phase genes in 40 cells raises their dataset means, so their
    // mean-matched controls sit slightly higher than the quiescent baseline
    // and the quiescent scores land a little below 0.)
    let quiescent_cycling = (40..120).filter(|&c| cc.cycling[c]).count();
    assert!(
        quiescent_cycling <= 80 / 2,
        "{quiescent_cycling} of 80 quiescent cells called cycling"
    );
    // Boosted cells score far above every quiescent cell.
    let max_quiescent_s = (40..120).map(|c| cc.s_score[c]).fold(f32::MIN, f32::max);
    let min_boosted_s = (0..20).map(|c| cc.s_score[c]).fold(f32::MAX, f32::min);
    assert!(
        min_boosted_s > max_quiescent_s + 0.5,
        "{min_boosted_s} vs {max_quiescent_s}"
    );
    for c in 40..120 {
        assert!(cc.s_score[c].abs() < 0.8, "{}", cc.s_score[c]);
        assert!(cc.g2m_score[c].abs() < 0.8, "{}", cc.g2m_score[c]);
    }
}

#[test]
fn undefined_without_cell_cycle_genes() {
    let m = Dense {
        genes: (0..30).map(|i| format!("G{i}")).collect(),
        cells: vec!["a".to_string(), "b".to_string()],
        counts: vec![vec![1, 2]; 30],
    };
    let cc = compute(&m, None);
    assert_eq!(cc.s_genes_mapped, 0);
    assert!(cc.s_score.iter().all(|v| v.is_nan()));
    assert!(cc.phase.iter().all(|p| *p == CellCyclePhase::Unknown));
    assert!(!cc.cycling.iter().any(|c| *c));
    assert_eq!(CellCyclePhase::Unknown.as_str(), "");
}

#[test]
fn legacy_symbols_are_recognised() {
    // MLF1IP (legacy) stands in for CENPU.
    let mut m = dataset(60);
    let idx = m.genes.iter().position(|g| g == "CENPU").unwrap();
    m.genes[idx] = "MLF1IP".to_string();
    let cc = compute(&m, None);
    assert_eq!(cc.s_genes_mapped, S_GENES.len());
}

#[test]
fn top2a_is_no_longer_a_conflict_risk_gene() {
    assert!(!CONFLICT_RISK_PANEL.contains(&"TOP2A"));
    assert!(G2M_GENES.contains(&"TOP2A"));
}

#[test]
fn species_is_inferred_from_symbol_casing() {
    use kira_spliceqc::genesets::aliases::detect_species;
    let human = dataset(10);
    assert_eq!(detect_species(&human), "human");
    let mut mouse = dataset(10);
    for g in mouse.genes.iter_mut() {
        let lower = g.to_ascii_lowercase();
        let mut chars = lower.chars();
        *g = chars.next().unwrap().to_ascii_uppercase().to_string() + chars.as_str();
    }
    assert_eq!(detect_species(&mouse), "mouse");
}
