use std::collections::HashMap;

use crate::expression::ExpressionMatrix;
use crate::model::splicing_instability::{
    PanelCoverage, RLOOP_RISK_HIGH_THRESHOLD, SPLICE_OVERLOAD_HIGH_THRESHOLD,
    SPLICING_INSTABILITY_HIGH_THRESHOLD, SplicingInstabilityGlobalStats,
    SplicingInstabilityMetrics, SplicingInstabilityMissingness, SplicingInstabilityRobustRef,
    SplicingInstabilityZReference,
};

use self::aggregate::aggregate_cluster_stats;
use self::panels::{
    CONFLICT_RISK_PANEL, MIN_GENES_PER_PANEL_CELL, NMD_PANEL, RLOOP_RESOLUTION_PANEL,
    SPLICEOSOME_PANEL, SPLICEQC_INSTABILITY_PANEL_V1, SPLICING_RBP_PANEL,
};
use self::scores::{panel_trimmed_mean, percentile, robust_z};

pub mod aggregate;
pub mod junction;
pub mod panels;
pub mod scores;

#[derive(Debug, Clone)]
struct ResolvedPanel {
    name: &'static str,
    genes_defined: usize,
    gene_indices: Vec<u32>,
}

pub fn compute(matrix: &dyn ExpressionMatrix) -> SplicingInstabilityMetrics {
    let n_cells = matrix.n_cells();
    let mut symbol_to_idx: HashMap<String, u32> = HashMap::with_capacity(matrix.n_genes());
    for gene_idx in 0..matrix.n_genes() {
        let symbol = matrix.gene_symbol(gene_idx).to_ascii_uppercase();
        symbol_to_idx.entry(symbol).or_insert(gene_idx as u32);
    }

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

    let mut splice_core = vec![f32::NAN; n_cells];
    let mut rbp_core = vec![f32::NAN; n_cells];
    let mut rloop_resolve_core = vec![f32::NAN; n_cells];
    let mut conflict_risk_core = vec![f32::NAN; n_cells];
    let mut nmd_core = vec![f32::NAN; n_cells];

    let mut scratch = Vec::with_capacity(64);

    for cell in 0..n_cells {
        splice_core[cell] = panel_trimmed_mean(
            matrix,
            &splice_panel.gene_indices,
            cell,
            MIN_GENES_PER_PANEL_CELL,
            &mut scratch,
        );
        rbp_core[cell] = panel_trimmed_mean(
            matrix,
            &rbp_panel.gene_indices,
            cell,
            MIN_GENES_PER_PANEL_CELL,
            &mut scratch,
        );
        rloop_resolve_core[cell] = panel_trimmed_mean(
            matrix,
            &rloop_panel.gene_indices,
            cell,
            MIN_GENES_PER_PANEL_CELL,
            &mut scratch,
        );

        if conflict_panel_enabled {
            conflict_risk_core[cell] = panel_trimmed_mean(
                matrix,
                &conflict_panel.gene_indices,
                cell,
                MIN_GENES_PER_PANEL_CELL,
                &mut scratch,
            );
        }

        if nmd_panel_enabled {
            nmd_core[cell] = panel_trimmed_mean(
                matrix,
                &nmd_panel.gene_indices,
                cell,
                MIN_GENES_PER_PANEL_CELL,
                &mut scratch,
            );
        }
    }

    let (z_splice, ref_splice) = robust_z(&splice_core);
    let (z_rbp, ref_rbp) = robust_z(&rbp_core);
    let (z_rloop, ref_rloop) = robust_z(&rloop_resolve_core);
    let (z_conflict, ref_conflict) = if conflict_panel_enabled {
        let (z, r) = robust_z(&conflict_risk_core);
        (z, Some(r))
    } else {
        (vec![0.0; n_cells], None)
    };
    let (z_nmd, ref_nmd) = if nmd_panel_enabled {
        let (z, r) = robust_z(&nmd_core);
        (z, Some(r))
    } else {
        (vec![0.0; n_cells], None)
    };

    let mut sos = vec![f32::NAN; n_cells];
    let mut rlr = vec![f32::NAN; n_cells];
    let mut sii = vec![f32::NAN; n_cells];

    let mut splice_overload_high = vec![false; n_cells];
    let mut rloop_risk_high = vec![false; n_cells];
    let mut splicing_instability_high = vec![false; n_cells];
    let mut genome_instability_splicing_flag = vec![false; n_cells];

    for cell in 0..n_cells {
        let zs = z_splice[cell];
        let zr = z_rbp[cell];
        if zs.is_finite() && zr.is_finite() {
            sos[cell] = 0.65 * zs + 0.35 * zr;
        }

        let zres = z_rloop[cell];
        rlr[cell] = if !zres.is_finite() {
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

        sii[cell] = if !sos[cell].is_finite() {
            f32::NAN
        } else if nmd_panel_enabled {
            let znmd = z_nmd[cell];
            if znmd.is_finite() {
                0.6 * relu(sos[cell]) + 0.4 * relu(-znmd)
            } else {
                f32::NAN
            }
        } else {
            relu(sos[cell])
        };

        splice_overload_high[cell] =
            sos[cell].is_finite() && sos[cell] >= SPLICE_OVERLOAD_HIGH_THRESHOLD;
        rloop_risk_high[cell] = rlr[cell].is_finite() && rlr[cell] >= RLOOP_RISK_HIGH_THRESHOLD;
        splicing_instability_high[cell] =
            sii[cell].is_finite() && sii[cell] >= SPLICING_INSTABILITY_HIGH_THRESHOLD;
        genome_instability_splicing_flag[cell] =
            splice_overload_high[cell] && rloop_risk_high[cell];
    }

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

fn resolve_panel(
    name: &'static str,
    genes: &[&str],
    symbol_to_idx: &HashMap<String, u32>,
) -> ResolvedPanel {
    let mut gene_indices = Vec::with_capacity(genes.len());
    for symbol in genes {
        if let Some(gene_idx) = symbol_to_idx.get(*symbol) {
            gene_indices.push(*gene_idx);
        }
    }
    ResolvedPanel {
        name,
        genes_defined: genes.len(),
        gene_indices,
    }
}

fn relu(value: f32) -> f32 {
    if value > 0.0 { value } else { 0.0 }
}

fn count_nan(values: &[f32]) -> usize {
    values.iter().filter(|v| !v.is_finite()).count()
}
