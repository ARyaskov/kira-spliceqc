//! Cell-cycle scoring (confounder guard).
//!
//! S-phase and G2/M scores follow Tirosh et al. 2016 Science (gene lists as
//! updated in Seurat's `cc.genes.updated.2019`): the mean `log1p(cp10k)` of
//! the phase genes minus the mean of control genes with matching expression
//! (the same `ControlPool` used for the splicing panels, with the cell-cycle
//! genes themselves excluded from the pool). Phase assignment follows Seurat
//! `CellCycleScoring`: S when the S score is the larger and positive, G2M
//! when the G2/M score is, else G1. Note that on a non-cycling tissue this
//! rule over-calls S/G2M (any positive noise counts), so the phase is a
//! confounder annotation, not a QC verdict.
//!
//! Splicing-factor and R-loop panels are enriched for genes that rise in
//! S/G2M (spliceosome, TOP2A, BRCA1/2); the `cycling` flag tells the reader
//! that a cell's splicing signatures may reflect proliferation.

use rayon::prelude::*;
use tracing::{info, warn};

use crate::expression::ExpressionMatrix;
use crate::genesets::aliases::{resolve_symbol, symbol_index};
use crate::genesets::controls::ControlPool;
use crate::model::cell_cycle::{CellCycleMetrics, CellCyclePhase};

/// Tirosh et al. 2016 S-phase genes (Seurat cc.genes.updated.2019$s.genes).
pub const S_GENES: &[&str] = &[
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM7", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7",
    "DTL", "PRIM1", "UHRF1", "CENPU", "HELLS", "RFC2", "POLR1B", "NASP", "RAD51AP1", "GMNN",
    "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45",
    "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B",
    "MRPL36", "E2F8",
];

/// Tirosh et al. 2016 G2/M genes (Seurat cc.genes.updated.2019$g2m.genes).
pub const G2M_GENES: &[&str] = &[
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2",
    "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "PIMREG", "SMC4", "CCNB2", "CKAP2L", "CKAP2",
    "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "JPT1",
    "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2",
    "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3",
    "GAS2L3", "CBX5", "CENPA",
];

/// Minimum mapped genes per list for the scores to be defined.
pub const MIN_GENES_PER_LIST: usize = 10;

/// Gene ids of every S / G2M gene present in `matrix` (for control-pool exclusion).
pub fn cell_cycle_gene_ids(matrix: &dyn ExpressionMatrix) -> Vec<u32> {
    let (s, g2m) = resolve(matrix);
    s.into_iter().chain(g2m).collect()
}

fn resolve(matrix: &dyn ExpressionMatrix) -> (Vec<u32>, Vec<u32>) {
    // Legacy symbols (MLF1IP, RPA2, FAM64A, HN1) resolve through the shared alias table.
    let index = symbol_index(matrix);
    let lookup = |symbols: &[&str]| -> Vec<u32> {
        let mut ids: Vec<u32> = symbols
            .iter()
            .filter_map(|s| resolve_symbol(&index, s).map(|(id, _)| id))
            .collect();
        ids.sort_unstable();
        ids.dedup();
        ids
    };
    (lookup(S_GENES), lookup(G2M_GENES))
}

pub fn compute(matrix: &dyn ExpressionMatrix, controls: Option<&ControlPool>) -> CellCycleMetrics {
    let n_cells = matrix.n_cells();
    let (s_ids, g2m_ids) = resolve(matrix);
    let defined = s_ids.len() >= MIN_GENES_PER_LIST && g2m_ids.len() >= MIN_GENES_PER_LIST;
    if !defined {
        warn!(
            s_genes_mapped = s_ids.len(),
            g2m_genes_mapped = g2m_ids.len(),
            min = MIN_GENES_PER_LIST,
            "too few cell-cycle genes mapped; phase scores undefined"
        );
        return CellCycleMetrics {
            s_genes_mapped: s_ids.len(),
            g2m_genes_mapped: g2m_ids.len(),
            s_score: vec![f32::NAN; n_cells],
            g2m_score: vec![f32::NAN; n_cells],
            phase: vec![CellCyclePhase::Unknown; n_cells],
            cycling: vec![false; n_cells],
        };
    }

    let s_ctrl = controls.map(|p| p.controls_for(&s_ids)).unwrap_or_default();
    let g2m_ctrl = controls.map(|p| p.controls_for(&g2m_ids)).unwrap_or_default();
    let score = |panel: &[u32], ctrl: &[u32], cell: usize| -> f32 {
        let scale = 1e4_f32 / matrix.libsize(cell).max(1) as f32;
        let mean = matrix.panel_ln1p_scaled_sum(panel, cell, scale) / panel.len() as f32;
        if ctrl.is_empty() {
            mean
        } else {
            mean - matrix.panel_ln1p_scaled_sum(ctrl, cell, scale) / ctrl.len() as f32
        }
    };

    let rows: Vec<(f32, f32)> = (0..n_cells)
        .into_par_iter()
        .map(|cell| (score(&s_ids, &s_ctrl, cell), score(&g2m_ids, &g2m_ctrl, cell)))
        .collect();

    let mut s_score = Vec::with_capacity(n_cells);
    let mut g2m_score = Vec::with_capacity(n_cells);
    let mut phase = Vec::with_capacity(n_cells);
    let mut cycling = Vec::with_capacity(n_cells);
    for (s, g) in rows {
        let ph = if !s.is_finite() || !g.is_finite() {
            CellCyclePhase::Unknown
        } else if s > g && s > 0.0 {
            CellCyclePhase::S
        } else if g > s && g > 0.0 {
            CellCyclePhase::G2M
        } else {
            CellCyclePhase::G1
        };
        cycling.push(matches!(ph, CellCyclePhase::S | CellCyclePhase::G2M));
        phase.push(ph);
        s_score.push(s);
        g2m_score.push(g);
    }
    info!(
        s_genes_mapped = s_ids.len(),
        g2m_genes_mapped = g2m_ids.len(),
        cycling = cycling.iter().filter(|c| **c).count(),
        "cell-cycle scores computed"
    );
    CellCycleMetrics {
        s_genes_mapped: s_ids.len(),
        g2m_genes_mapped: g2m_ids.len(),
        s_score,
        g2m_score,
        phase,
        cycling,
    }
}
