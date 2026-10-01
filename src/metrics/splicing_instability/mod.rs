use ahash::AHashMap;
use rayon::prelude::*;

use crate::expression::ExpressionMatrix;
use crate::genesets::aliases::{resolve_symbol, symbol_index};
use crate::genesets::controls::ControlPool;
use crate::model::splicing_instability::{
    PanelCoverage, SplicingInstabilityGlobalStats, SplicingInstabilityMetrics,
    SplicingInstabilityMissingness, SplicingInstabilityRobustRef, SplicingInstabilityZReference,
};

use self::aggregate::aggregate_cluster_stats;
use self::panels::{
    CONFLICT_RISK_PANEL, MIN_GENES_PER_PANEL_CELL, NMD_PANEL, RLOOP_RESOLUTION_PANEL,
    SPLICEOSOME_PANEL, SPLICEQC_INSTABILITY_PANEL_V1, SPLICING_RBP_PANEL,
};
use self::scores::{panel_trimmed_mean, percentile};
use crate::reference::{Strata, flag_outliers, robust_z_by_stratum, robust_z_by_stratum_and_depth};
use crate::stats::robust::{RobustRef, robust_z_logged};

pub mod aggregate;
pub mod junction;
pub mod panels;
pub mod scores;

#[derive(Debug, Clone)]
struct ResolvedPanel {
    name: &'static str,
    genes_defined: usize,
    /// Sorted ascending — required for sparse-merge gather.
    gene_indices: Vec<u32>,
}

/// Panel cores are trimmed means of `log1p(cp10k)`; with a `ControlPool` the
/// mean of each panel's control set is subtracted (depth correction).
pub fn compute(
    matrix: &dyn ExpressionMatrix,
    controls: Option<&ControlPool>,
    strata: &Strata,
) -> SplicingInstabilityMetrics {
    let n_cells = matrix.n_cells();

    // Case-insensitive index with legacy aliases, first occurrence wins
    // (same resolution as the geneset loader).
    let symbol_to_idx = symbol_index(matrix);

    let splice_panel = resolve_panel("spliceosome_panel", SPLICEOSOME_PANEL, &symbol_to_idx);
    let rbp_panel = resolve_panel("splicing_rbp_panel", SPLICING_RBP_PANEL, &symbol_to_idx);
    let rloop_panel = resolve_panel(
        "rloop_resolution_panel",
        RLOOP_RESOLUTION_PANEL,
        &symbol_to_idx,
    );
    let conflict_panel = resolve_panel("conflict_risk_panel", CONFLICT_RISK_PANEL, &symbol_to_idx);
    let nmd_panel = resolve_panel("nmd_panel", NMD_PANEL, &symbol_to_idx);

    let conflict_panel_enabled = conflict_panel.gene_indices.len() >= MIN_GENES_PER_PANEL_CELL;
    let nmd_panel_enabled = nmd_panel.gene_indices.len() >= MIN_GENES_PER_PANEL_CELL;

    let control_sets: Vec<Vec<u32>> = [
        &splice_panel,
        &rbp_panel,
        &rloop_panel,
        &conflict_panel,
        &nmd_panel,
    ]
    .iter()
    .map(|p| {
        controls
            .map(|c| c.controls_for(&p.gene_indices))
            .unwrap_or_default()
    })
    .collect();
    let libsize_scale: Vec<f32> = (0..n_cells)
        .map(|cell| 1e4_f32 / matrix.libsize(cell).max(1) as f32)
        .collect();
    let background = |k: usize, cell: usize| -> f32 {
        let set = &control_sets[k];
        if set.is_empty() {
            0.0
        } else {
            matrix.panel_ln1p_scaled_sum(set, cell, libsize_scale[cell]) / set.len() as f32
        }
    };

    // One parallel per-cell pass over all enabled panels (each gets its own scratch buf).
    let panel_rows: Vec<[f32; 5]> = (0..n_cells)
        .into_par_iter()
        .map_init(
            || Vec::<f32>::with_capacity(64),
            |scratch, cell| {
                let splice = panel_trimmed_mean(
                    matrix,
                    &splice_panel.gene_indices,
                    cell,
                    MIN_GENES_PER_PANEL_CELL,
                    scratch,
                ) - background(0, cell);
                let rbp = panel_trimmed_mean(
                    matrix,
                    &rbp_panel.gene_indices,
                    cell,
                    MIN_GENES_PER_PANEL_CELL,
                    scratch,
                ) - background(1, cell);
                let rloop = panel_trimmed_mean(
                    matrix,
                    &rloop_panel.gene_indices,
                    cell,
                    MIN_GENES_PER_PANEL_CELL,
                    scratch,
                ) - background(2, cell);
                let conflict = if conflict_panel_enabled {
                    panel_trimmed_mean(
                        matrix,
                        &conflict_panel.gene_indices,
                        cell,
                        MIN_GENES_PER_PANEL_CELL,
                        scratch,
                    ) - background(3, cell)
                } else {
                    f32::NAN
                };
                let nmd = if nmd_panel_enabled {
                    panel_trimmed_mean(
                        matrix,
                        &nmd_panel.gene_indices,
                        cell,
                        MIN_GENES_PER_PANEL_CELL,
                        scratch,
                    ) - background(4, cell)
                } else {
                    f32::NAN
                };
                [splice, rbp, rloop, conflict, nmd]
            },
        )
        .collect();

    let mut splice_core = Vec::with_capacity(n_cells);
    let mut rbp_core = Vec::with_capacity(n_cells);
    let mut rloop_resolve_core = Vec::with_capacity(n_cells);
    let mut conflict_risk_core = Vec::with_capacity(n_cells);
    let mut nmd_core = Vec::with_capacity(n_cells);
    for row in panel_rows {
        splice_core.push(row[0]);
        rbp_core.push(row[1]);
        rloop_resolve_core.push(row[2]);
        conflict_risk_core.push(row[3]);
        nmd_core.push(row[4]);
    }

    // z-scores are standardized within stratum and library-size bin; the
    // reported `z_reference` medians/MADs are the dataset-wide values kept
    // for provenance and cross-sample comparison.
    let libsize: Vec<u64> = (0..n_cells).map(|c| matrix.libsize(c)).collect();
    let standardize = |values: &[f32], name: &str| -> (Vec<f32>, RobustRef) {
        let (_, r) = robust_z_logged(values, name);
        let (z, _) = robust_z_by_stratum_and_depth(values, strata, &libsize);
        (z, r)
    };
    let (z_splice, ref_splice) = standardize(&splice_core, "splice_core");
    let (z_rbp, ref_rbp) = standardize(&rbp_core, "rbp_core");
    let (z_rloop, ref_rloop) = standardize(&rloop_resolve_core, "rloop_resolve_core");
    let (z_conflict, ref_conflict) = if conflict_panel_enabled {
        let (z, r) = standardize(&conflict_risk_core, "conflict_risk_core");
        (z, Some(r))
    } else {
        (vec![0.0; n_cells], None)
    };
    let (z_nmd, ref_nmd) = if nmd_panel_enabled {
        let (z, r) = standardize(&nmd_core, "nmd_core");
        (z, Some(r))
    } else {
        (vec![0.0; n_cells], None)
    };

    // Composite scores — independent per cell.
    // Composites keep their documented (relu-based) definitions; the flags
    // come from signed versions standardized within the stratum, so they are
    // calibrated like every other flag instead of using fixed cut-offs.
    let composite: Vec<(f32, f32, f32, f32, f32, f32)> = (0..n_cells)
        .into_par_iter()
        .map(|cell| {
            let zs = z_splice[cell];
            let zr = z_rbp[cell];
            let sos = if zs.is_finite() && zr.is_finite() {
                0.65 * zs + 0.35 * zr
            } else {
                f32::NAN
            };

            let zres = z_rloop[cell];
            let rlr = if !zres.is_finite() {
                f32::NAN
            } else if conflict_panel_enabled {
                let zrisk = z_conflict[cell];
                if zrisk.is_finite() {
                    0.7 * relu(-zres) + 0.3 * relu(zrisk)
                } else {
                    f32::NAN
                }
            } else {
                relu(-zres)
            };

            let sii = if !sos.is_finite() {
                f32::NAN
            } else if nmd_panel_enabled {
                let znmd = z_nmd[cell];
                if znmd.is_finite() {
                    0.6 * relu(sos) + 0.4 * relu(-znmd)
                } else {
                    f32::NAN
                }
            } else {
                relu(sos)
            };

            // Signed counterparts (no relu) for calibrated flags.
            let rlr_signed = if !zres.is_finite() {
                f32::NAN
            } else if conflict_panel_enabled {
                let zrisk = z_conflict[cell];
                if zrisk.is_finite() {
                    0.7 * -zres + 0.3 * zrisk
                } else {
                    f32::NAN
                }
            } else {
                -zres
            };
            let sii_signed = if !sos.is_finite() {
                f32::NAN
            } else if nmd_panel_enabled {
                let znmd = z_nmd[cell];
                if znmd.is_finite() {
                    0.6 * sos + 0.4 * -znmd
                } else {
                    f32::NAN
                }
            } else {
                sos
            };

            (sos, rlr, sii, sos, rlr_signed, sii_signed)
        })
        .collect();

    let mut sos = Vec::with_capacity(n_cells);
    let mut rlr = Vec::with_capacity(n_cells);
    let mut sii = Vec::with_capacity(n_cells);
    let mut sos_signed = Vec::with_capacity(n_cells);
    let mut rlr_signed = Vec::with_capacity(n_cells);
    let mut sii_signed = Vec::with_capacity(n_cells);
    for (s, r, i, ss, rs, is) in composite {
        sos.push(s);
        rlr.push(r);
        sii.push(i);
        sos_signed.push(ss);
        rlr_signed.push(rs);
        sii_signed.push(is);
    }
    let (sos_dev, _) = robust_z_by_stratum(&sos_signed, strata);
    let (rlr_dev, _) = robust_z_by_stratum(&rlr_signed, strata);
    let (sii_dev, _) = robust_z_by_stratum(&sii_signed, strata);
    let splice_overload_high = flag_outliers(&sos_dev, strata, 1.0);
    let rloop_risk_high = flag_outliers(&rlr_dev, strata, 1.0);
    let splicing_instability_high = flag_outliers(&sii_dev, strata, 1.0);
    let genome_instability_splicing_flag: Vec<bool> = splice_overload_high
        .iter()
        .zip(&rloop_risk_high)
        .map(|(a, b)| *a && *b)
        .collect();

    let global_stats = SplicingInstabilityGlobalStats {
        sos_p50: percentile(&sos, 0.5),
        sos_p90: percentile(&sos, 0.9),
        rlr_p50: percentile(&rlr, 0.5),
        rlr_p90: percentile(&rlr, 0.9),
        sii_p50: percentile(&sii, 0.5),
        sii_p90: percentile(&sii, 0.9),
    };

    let z_reference = SplicingInstabilityZReference {
        splice_core: SplicingInstabilityRobustRef {
            median: ref_splice.median,
            mad: ref_splice.mad,
        },
        rbp_core: SplicingInstabilityRobustRef {
            median: ref_rbp.median,
            mad: ref_rbp.mad,
        },
        rloop_resolve_core: SplicingInstabilityRobustRef {
            median: ref_rloop.median,
            mad: ref_rloop.mad,
        },
        conflict_risk_core: ref_conflict.map(|r| SplicingInstabilityRobustRef {
            median: r.median,
            mad: r.mad,
        }),
        nmd_core: ref_nmd.map(|r| SplicingInstabilityRobustRef {
            median: r.median,
            mad: r.mad,
        }),
    };

    let missingness = SplicingInstabilityMissingness {
        splice_core_nan_cells: count_nan(&splice_core),
        rbp_core_nan_cells: count_nan(&rbp_core),
        rloop_resolve_core_nan_cells: count_nan(&rloop_resolve_core),
        conflict_risk_core_nan_cells: count_nan(&conflict_risk_core),
        nmd_core_nan_cells: count_nan(&nmd_core),
        sos_nan_cells: count_nan(&sos),
        rlr_nan_cells: count_nan(&rlr),
        sii_nan_cells: count_nan(&sii),
        panel_coverage: vec![
            PanelCoverage {
                panel_name: splice_panel.name,
                genes_defined: splice_panel.genes_defined,
                genes_mapped: splice_panel.gene_indices.len(),
            },
            PanelCoverage {
                panel_name: rbp_panel.name,
                genes_defined: rbp_panel.genes_defined,
                genes_mapped: rbp_panel.gene_indices.len(),
            },
            PanelCoverage {
                panel_name: rloop_panel.name,
                genes_defined: rloop_panel.genes_defined,
                genes_mapped: rloop_panel.gene_indices.len(),
            },
            PanelCoverage {
                panel_name: conflict_panel.name,
                genes_defined: conflict_panel.genes_defined,
                genes_mapped: conflict_panel.gene_indices.len(),
            },
            PanelCoverage {
                panel_name: nmd_panel.name,
                genes_defined: nmd_panel.genes_defined,
                genes_mapped: nmd_panel.gene_indices.len(),
            },
        ],
    };

    let cluster_stats = aggregate_cluster_stats(
        None,
        &sos,
        &rlr,
        &sii,
        &splice_overload_high,
        &rloop_risk_high,
        &splicing_instability_high,
        &genome_instability_splicing_flag,
    );

    SplicingInstabilityMetrics {
        panel_version: SPLICEQC_INSTABILITY_PANEL_V1,
        min_genes: MIN_GENES_PER_PANEL_CELL,
        conflict_panel_enabled,
        nmd_panel_enabled,
        splice_core,
        rbp_core,
        rloop_resolve_core,
        conflict_risk_core,
        nmd_core,
        sos,
        rlr,
        sii,
        sos_dev,
        rlr_dev,
        sii_dev,
        splice_overload_high,
        rloop_risk_high,
        splicing_instability_high,
        genome_instability_splicing_flag,
        z_reference,
        global_stats,
        cluster_stats,
        missingness,
    }
}

/// Gene ids (in `matrix`) of every stage-15 panel gene, for control-pool exclusion.
pub fn panel_gene_ids(matrix: &dyn ExpressionMatrix) -> Vec<u32> {
    let symbol_to_idx = symbol_index(matrix);
    let mut ids = Vec::new();
    for panel in [
        SPLICEOSOME_PANEL,
        SPLICING_RBP_PANEL,
        RLOOP_RESOLUTION_PANEL,
        CONFLICT_RISK_PANEL,
        NMD_PANEL,
    ] {
        ids.extend(resolve_panel("", panel, &symbol_to_idx).gene_indices);
    }
    ids
}

fn resolve_panel(
    name: &'static str,
    genes: &[&str],
    symbol_to_idx: &AHashMap<String, u32>,
) -> ResolvedPanel {
    let mut gene_indices = Vec::with_capacity(genes.len());
    for symbol in genes {
        if let Some((gene_idx, _)) = resolve_symbol(symbol_to_idx, symbol) {
            gene_indices.push(gene_idx);
        }
    }
    // Sort for sparse-merge gather. Dedup defensive (panel literals shouldn't repeat).
    gene_indices.sort_unstable();
    gene_indices.dedup();
    ResolvedPanel {
        name,
        genes_defined: genes.len(),
        gene_indices,
    }
}

#[inline]
fn relu(value: f32) -> f32 {
    if value > 0.0 { value } else { 0.0 }
}

fn count_nan(values: &[f32]) -> usize {
    values.iter().filter(|v| !v.is_finite()).count()
}
