use std::path::Path;

use serde::Serialize;

use crate::input::error::InputError;
use crate::model::assembly_phase::AssemblyPhaseImbalanceMetrics;
use crate::model::collapse::{SpliceosomeCollapseMetrics, SpliceosomeCollapseStatus};
use crate::model::coupling::CouplingStressMetrics;
use crate::model::cryptic_risk::CrypticSplicingRiskMetrics;
use crate::model::exon_intron_bias::ExonIntronDefinitionMetrics;
use crate::model::imbalance::SpliceosomeImbalanceMetrics;
use crate::model::isoform_dispersion::IsoformDispersionMetrics;
use crate::model::missplicing::MissplicingMetrics;
use crate::model::sis::{SpliceIntegrityClass, SpliceIntegrityMetrics};
use crate::model::splicing_instability::{
    RLOOP_RISK_HIGH_THRESHOLD, SPLICE_OVERLOAD_HIGH_THRESHOLD, SPLICING_INSTABILITY_HIGH_THRESHOLD,
    SplicingInstabilityMetrics,
};
use crate::model::splicing_noise::SplicingNoiseMetrics;
use crate::model::timecourse::{SplicingTrajectoryClass, TimecourseSplicingMetrics};
use crate::model::unspliced::UnsplicedMetrics;

/// Bumped to 2.0 in v0.3: expression-signature keys carry the `_expr`
/// suffix (`regulator_expr`, `missplicing_expr`, `imbalance_expr`,
/// `coupling_expr`, `exon_intron_bias_expr`, `assembly_phase_expr`, and the
/// panel cores inside `splicing_instability`).
pub const JSON_SCHEMA_VERSION: &str = "2.0";

#[derive(Serialize)]
struct JsonOutput<'a> {
    schema_version: &'static str,
    tool: &'static str,
    mode: &'static str,
    n_cells: usize,
    /// Whether experimental composite signatures (sis/class, SOS/RLR/SII,
    /// flags, cryptic risk, collapse) are included.
    experimental_signatures: bool,
    cells: Vec<JsonCell<'a>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    sis: Option<JsonSisStage<'a>>,
    #[serde(rename = "regulator_expr")]
    isoform: JsonIsoformStage,
    #[serde(rename = "missplicing_expr")]
    missplicing: JsonMissplicingStage,
    #[serde(rename = "imbalance_expr")]
    imbalance: JsonImbalanceStage,
    splicing_noise: Option<JsonSplicingNoise>,
    cryptic_risk: Option<JsonCrypticRisk>,
    collapse: Option<JsonCollapse>,
    timecourse: Option<JsonTimecourse>,
    splicing_instability: JsonSplicingInstabilityStage,
    cell_cycle_guardrail: Option<JsonCellCycleGuardrail>,
    /// Input levels available: `["L0"]` or `["L0", "L1"]`.
    input_levels: Vec<&'static str>,
    /// Tier A unspliced metrics (input level L1); absent without layers.
    #[serde(skip_serializing_if = "Option::is_none")]
    unspliced: Option<JsonUnsplicedStage>,
}

#[derive(Serialize)]
struct JsonCellUnspliced {
    spliced_umis: u64,
    unspliced_umis: u64,
    ambiguous_umis: u64,
    unspliced_fraction: Option<f32>,
    unspliced_fraction_ci_low: Option<f32>,
    unspliced_fraction_ci_high: Option<f32>,
}

#[derive(Serialize)]
struct JsonUnsplicedStage {
    source: String,
    min_layer_umis: u64,
    has_ambiguous: bool,
    cells_without_layers: usize,
    undefined_cells: usize,
    spliced_umis: Vec<u64>,
    unspliced_umis: Vec<u64>,
    ambiguous_umis: Vec<u64>,
    unspliced_fraction: Vec<Option<f32>>,
    unspliced_fraction_ci_low: Vec<Option<f32>>,
    unspliced_fraction_ci_high: Vec<Option<f32>>,
}

#[derive(Serialize)]
struct JsonCell<'a> {
    cell_id: usize,
    cell_name: &'a str,
    #[serde(skip_serializing_if = "Option::is_none")]
    sis: Option<f32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    class: Option<&'a str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    penalties: Option<JsonPenalties>,
    #[serde(rename = "regulator_expr")]
    isoform: JsonIsoform,
    #[serde(rename = "missplicing_expr")]
    missplicing: JsonMissplicing,
    #[serde(rename = "imbalance_expr")]
    imbalance: JsonImbalance,
    #[serde(rename = "coupling_expr")]
    coupling: JsonCoupling,
    #[serde(rename = "exon_intron_bias_expr")]
    exon_intron_bias: JsonExonIntronBias,
    #[serde(rename = "assembly_phase_expr")]
    assembly_phase: JsonAssemblyPhase,
    splicing_instability: JsonCellSplicingInstability,
    #[serde(skip_serializing_if = "Option::is_none")]
    unspliced: Option<JsonCellUnspliced>,
}

#[derive(Serialize)]
struct JsonPenalties {
    missplicing: Option<f32>,
    imbalance: Option<f32>,
    entropy_z: Option<f32>,
    entropy_abs: Option<f32>,
}

#[derive(Serialize)]
struct JsonIsoform {
    entropy: Option<f32>,
    dispersion: Option<f32>,
    z_entropy: Option<f32>,
}

#[derive(Serialize)]
struct JsonMissplicing {
    core: Option<f32>,
    u12: Option<f32>,
    nmd: Option<f32>,
    sr_hnrnp: Option<f32>,
    burden: Option<f32>,
}

#[derive(Serialize)]
struct JsonImbalance {
    z_u1: Option<f32>,
    z_u2: Option<f32>,
    z_sf3b: Option<f32>,
    z_srsf: Option<f32>,
    z_hnrnp: Option<f32>,
    z_u12: Option<f32>,
    z_nmd: Option<f32>,
    axis_sr_hnrnp: Option<f32>,
    axis_u2_u1: Option<f32>,
    axis_u12_major: Option<f32>,
    axis_nmd: Option<f32>,
    imbalance: Option<f32>,
}

#[derive(Serialize)]
struct JsonCoupling {
    coupling_stress: Option<f32>,
}

#[derive(Serialize)]
struct JsonExonIntronBias {
    exon_definition_bias: Option<f32>,
    z_srsf: Option<f32>,
    z_u2af: Option<f32>,
    z_hnrnp: Option<f32>,
}

#[derive(Serialize)]
struct JsonAssemblyPhase {
    z_ea: Option<f32>,
    z_b: Option<f32>,
    z_cat: Option<f32>,
    ea_imbalance: Option<f32>,
    b_imbalance: Option<f32>,
    cat_imbalance: Option<f32>,
}

#[derive(Serialize)]
struct JsonCellSplicingInstability {
    #[serde(rename = "spliceosome_core_expr")]
    splice_core: Option<f32>,
    #[serde(rename = "splicing_rbp_expr")]
    rbp_core: Option<f32>,
    #[serde(rename = "rloop_resolution_expr")]
    rloop_resolve_core: Option<f32>,
    #[serde(rename = "conflict_risk_expr")]
    conflict_risk_core: Option<f32>,
    #[serde(rename = "nmd_factor_expr")]
    nmd_core: Option<f32>,
    #[serde(flatten, skip_serializing_if = "Option::is_none")]
    composites: Option<JsonCellSplicingComposites>,
}

/// Experimental composite scores and flags (only with `experimental`).
#[derive(Serialize)]
struct JsonCellSplicingComposites {
    sos: Option<f32>,
    rlr: Option<f32>,
    sii: Option<f32>,
    splice_overload_high: bool,
    rloop_risk_high: bool,
    splicing_instability_high: bool,
    genome_instability_splicing_flag: bool,
}

#[derive(Serialize)]
struct JsonSisStage<'a> {
    sis: Vec<Option<f32>>,
    class: Vec<&'a str>,
    p_missplicing: Vec<Option<f32>>,
    p_imbalance: Vec<Option<f32>>,
    p_entropy_z: Vec<Option<f32>>,
    p_entropy_abs: Vec<Option<f32>>,
}

#[derive(Serialize)]
struct JsonIsoformStage {
    entropy: Vec<Option<f32>>,
    dispersion: Vec<Option<f32>>,
    z_entropy: Vec<Option<f32>>,
}

#[derive(Serialize)]
struct JsonMissplicingStage {
    b_core: Vec<Option<f32>>,
    b_u12: Vec<Option<f32>>,
    b_nmd: Vec<Option<f32>>,
    b_srhn: Vec<Option<f32>>,
    burden: Vec<Option<f32>>,
    burden_star: Vec<Option<f32>>,
}

#[derive(Serialize)]
struct JsonImbalanceStage {
    z_u1: Vec<Option<f32>>,
    z_u2: Vec<Option<f32>>,
    z_sf3b: Vec<Option<f32>>,
    z_srsf: Vec<Option<f32>>,
    z_hnrnp: Vec<Option<f32>>,
    z_u12: Vec<Option<f32>>,
    z_nmd: Vec<Option<f32>>,
    axis_sr_hnrnp: Vec<Option<f32>>,
    axis_u2_u1: Vec<Option<f32>>,
    axis_u12_major: Vec<Option<f32>>,
    axis_nmd: Vec<Option<f32>>,
    imbalance: Vec<Option<f32>>,
}

#[derive(Serialize)]
struct JsonGenesetNoise {
    geneset: String,
    noise: Option<f32>,
}

#[derive(Serialize)]
struct JsonSplicingNoise {
    noise_index: Option<f32>,
    per_geneset_noise: Vec<JsonGenesetNoise>,
}

#[derive(Serialize)]
struct JsonCrypticRisk {
    cryptic_risk: Vec<Option<f32>>,
    x_sr_hnrnp: Vec<Option<f32>>,
    x_entropy: Vec<Option<f32>>,
    x_nmd: Vec<Option<f32>>,
}

#[derive(Serialize)]
struct JsonCollapse {
    collapse_status: Vec<&'static str>,
    core_suppression: Vec<bool>,
    high_imbalance: Vec<bool>,
    low_sis: Vec<bool>,
}

#[derive(Serialize)]
struct JsonTimecourse {
    trajectory: &'static str,
    delta_sis: Vec<Option<f32>>,
    delta_entropy: Vec<Option<f32>>,
    delta_imbalance: Vec<Option<f32>>,
}

#[derive(Serialize)]
struct JsonCellCycleGuardrail {}

#[derive(Serialize)]
struct JsonPanelCoverage {
    panel_name: &'static str,
    genes_defined: usize,
    genes_mapped: usize,
}

#[derive(Serialize)]
struct JsonSplicingInstabilityMissingness {
    splice_core_nan_cells: usize,
    rbp_core_nan_cells: usize,
    rloop_resolve_core_nan_cells: usize,
    conflict_risk_core_nan_cells: usize,
    nmd_core_nan_cells: usize,
    sos_nan_cells: usize,
    rlr_nan_cells: usize,
    sii_nan_cells: usize,
    panel_coverage: Vec<JsonPanelCoverage>,
}

#[derive(Serialize)]
struct JsonRobustRef {
    median: Option<f32>,
    mad: Option<f32>,
}

#[derive(Serialize)]
struct JsonSplicingInstabilityZReference {
    #[serde(rename = "spliceosome_core_expr")]
    splice_core: JsonRobustRef,
    #[serde(rename = "splicing_rbp_expr")]
    rbp_core: JsonRobustRef,
    #[serde(rename = "rloop_resolution_expr")]
    rloop_resolve_core: JsonRobustRef,
    #[serde(rename = "conflict_risk_expr")]
    conflict_risk_core: Option<JsonRobustRef>,
    #[serde(rename = "nmd_factor_expr")]
    nmd_core: Option<JsonRobustRef>,
}

#[derive(Serialize)]
struct JsonSplicingInstabilityGlobalStats {
    sos_p50: Option<f32>,
    sos_p90: Option<f32>,
    rlr_p50: Option<f32>,
    rlr_p90: Option<f32>,
    sii_p50: Option<f32>,
    sii_p90: Option<f32>,
}

#[derive(Serialize)]
struct JsonSplicingInstabilityThresholds {
    splice_overload_high: f32,
    rloop_risk_high: f32,
    splicing_instability_high: f32,
}

#[derive(Serialize)]
struct JsonSplicingInstabilityStage {
    panel_version: &'static str,
    min_genes_per_panel_cell: usize,
    conflict_panel_enabled: bool,
    nmd_panel_enabled: bool,
    #[serde(rename = "spliceosome_core_expr")]
    splice_core: Vec<Option<f32>>,
    #[serde(rename = "splicing_rbp_expr")]
    rbp_core: Vec<Option<f32>>,
    #[serde(rename = "rloop_resolution_expr")]
    rloop_resolve_core: Vec<Option<f32>>,
    #[serde(rename = "conflict_risk_expr")]
    conflict_risk_core: Vec<Option<f32>>,
    #[serde(rename = "nmd_factor_expr")]
    nmd_core: Vec<Option<f32>>,
    z_reference: JsonSplicingInstabilityZReference,
    missingness: JsonSplicingInstabilityMissingness,
    #[serde(flatten, skip_serializing_if = "Option::is_none")]
    composites: Option<JsonSplicingInstabilityComposites>,
}

/// Experimental composite arrays, thresholds and their percentiles.
#[derive(Serialize)]
struct JsonSplicingInstabilityComposites {
    thresholds: JsonSplicingInstabilityThresholds,
    sos: Vec<Option<f32>>,
    rlr: Vec<Option<f32>>,
    sii: Vec<Option<f32>>,
    splice_overload_high: Vec<bool>,
    rloop_risk_high: Vec<bool>,
    splicing_instability_high: Vec<bool>,
    genome_instability_splicing_flag: Vec<bool>,
    global_stats: JsonSplicingInstabilityGlobalStats,
}

#[allow(clippy::too_many_arguments)]
pub fn write_json(
    path: &Path,
    cell_names: &[String],
    isoform: &IsoformDispersionMetrics,
    missplicing: &MissplicingMetrics,
    imbalance: &SpliceosomeImbalanceMetrics,
    sis: &SpliceIntegrityMetrics,
    coupling: Option<&CouplingStressMetrics>,
    exon_intron: Option<&ExonIntronDefinitionMetrics>,
    assembly: Option<&AssemblyPhaseImbalanceMetrics>,
    splicing_noise: Option<&SplicingNoiseMetrics>,
    cryptic_risk: Option<&CrypticSplicingRiskMetrics>,
    collapse: Option<&SpliceosomeCollapseMetrics>,
    timecourse: Option<&TimecourseSplicingMetrics>,
    splicing_instability: &SplicingInstabilityMetrics,
    unspliced: Option<&UnsplicedMetrics>,
    experimental: bool,
) -> Result<(), InputError> {
    let n_cells = cell_names.len();
    let mut cells = Vec::with_capacity(n_cells);

    // Composite stages are experimental: dropped from the output unless asked for.
    let cryptic_risk = if experimental { cryptic_risk } else { None };
    let collapse = if experimental { collapse } else { None };

    let coupling_default;
    let exon_intron_default;
    let assembly_default;

    let coupling_ref = match coupling {
        Some(value) => value,
        None => {
            coupling_default = CouplingStressMetrics {
                coupling_stress: vec![f32::NAN; n_cells],
            };
            &coupling_default
        }
    };

    let exon_intron_ref = match exon_intron {
        Some(value) => value,
        None => {
            exon_intron_default = ExonIntronDefinitionMetrics {
                exon_definition_bias: vec![f32::NAN; n_cells],
                z_srsf: vec![f32::NAN; n_cells],
                z_u2af: vec![f32::NAN; n_cells],
                z_hnrnp: vec![f32::NAN; n_cells],
            };
            &exon_intron_default
        }
    };

    let assembly_ref = match assembly {
        Some(value) => value,
        None => {
            assembly_default = AssemblyPhaseImbalanceMetrics {
                z_ea: vec![f32::NAN; n_cells],
                z_b: vec![f32::NAN; n_cells],
                z_cat: vec![f32::NAN; n_cells],
                ea_imbalance: vec![f32::NAN; n_cells],
                b_imbalance: vec![f32::NAN; n_cells],
                cat_imbalance: vec![f32::NAN; n_cells],
            };
            &assembly_default
        }
    };

    for (cell_id, cell_name) in cell_names.iter().enumerate() {
        let (sis_value, class, penalties) = if experimental {
            (
                opt_f32(sis.sis[cell_id]),
                Some(class_str(sis.class[cell_id])),
                Some(JsonPenalties {
                    missplicing: opt_f32(sis.p_missplicing[cell_id]),
                    imbalance: opt_f32(sis.p_imbalance[cell_id]),
                    entropy_z: opt_f32(sis.p_entropy_z[cell_id]),
                    entropy_abs: opt_f32(sis.p_entropy_abs[cell_id]),
                }),
            )
        } else {
            (None, None, None)
        };
        cells.push(JsonCell {
            cell_id,
            cell_name,
            sis: sis_value,
            class,
            penalties,
            isoform: JsonIsoform {
                entropy: opt_f32(isoform.entropy[cell_id]),
                dispersion: opt_f32(isoform.dispersion[cell_id]),
                z_entropy: opt_f32(isoform.z_entropy[cell_id]),
            },
            missplicing: JsonMissplicing {
                core: opt_f32(missplicing.b_core[cell_id]),
                u12: opt_f32(missplicing.b_u12[cell_id]),
                nmd: opt_f32(missplicing.b_nmd[cell_id]),
                sr_hnrnp: opt_f32(missplicing.b_srhn[cell_id]),
                burden: opt_f32(missplicing.burden[cell_id]),
            },
            imbalance: JsonImbalance {
                z_u1: opt_f32(imbalance.z_u1[cell_id]),
                z_u2: opt_f32(imbalance.z_u2[cell_id]),
                z_sf3b: opt_f32(imbalance.z_sf3b[cell_id]),
                z_srsf: opt_f32(imbalance.z_srsf[cell_id]),
                z_hnrnp: opt_f32(imbalance.z_hnrnp[cell_id]),
                z_u12: opt_f32(imbalance.z_u12[cell_id]),
                z_nmd: opt_f32(imbalance.z_nmd[cell_id]),
                axis_sr_hnrnp: opt_f32(imbalance.axis_sr_hnrnp[cell_id]),
                axis_u2_u1: opt_f32(imbalance.axis_u2_u1[cell_id]),
                axis_u12_major: opt_f32(imbalance.axis_u12_major[cell_id]),
                axis_nmd: opt_f32(imbalance.axis_nmd[cell_id]),
                imbalance: opt_f32(imbalance.imbalance[cell_id]),
            },
            coupling: JsonCoupling {
                coupling_stress: opt_f32(coupling_ref.coupling_stress[cell_id]),
            },
            exon_intron_bias: JsonExonIntronBias {
                exon_definition_bias: opt_f32(exon_intron_ref.exon_definition_bias[cell_id]),
                z_srsf: opt_f32(exon_intron_ref.z_srsf[cell_id]),
                z_u2af: opt_f32(exon_intron_ref.z_u2af[cell_id]),
                z_hnrnp: opt_f32(exon_intron_ref.z_hnrnp[cell_id]),
            },
            assembly_phase: JsonAssemblyPhase {
                z_ea: opt_f32(assembly_ref.z_ea[cell_id]),
                z_b: opt_f32(assembly_ref.z_b[cell_id]),
                z_cat: opt_f32(assembly_ref.z_cat[cell_id]),
                ea_imbalance: opt_f32(assembly_ref.ea_imbalance[cell_id]),
                b_imbalance: opt_f32(assembly_ref.b_imbalance[cell_id]),
                cat_imbalance: opt_f32(assembly_ref.cat_imbalance[cell_id]),
            },
            splicing_instability: JsonCellSplicingInstability {
                splice_core: opt_f32(splicing_instability.splice_core[cell_id]),
                rbp_core: opt_f32(splicing_instability.rbp_core[cell_id]),
                rloop_resolve_core: opt_f32(splicing_instability.rloop_resolve_core[cell_id]),
                conflict_risk_core: opt_f32(splicing_instability.conflict_risk_core[cell_id]),
                nmd_core: opt_f32(splicing_instability.nmd_core[cell_id]),
                composites: experimental.then(|| JsonCellSplicingComposites {
                    sos: opt_f32(splicing_instability.sos[cell_id]),
                    rlr: opt_f32(splicing_instability.rlr[cell_id]),
                    sii: opt_f32(splicing_instability.sii[cell_id]),
                    splice_overload_high: splicing_instability.splice_overload_high[cell_id],
                    rloop_risk_high: splicing_instability.rloop_risk_high[cell_id],
                    splicing_instability_high: splicing_instability.splicing_instability_high
                        [cell_id],
                    genome_instability_splicing_flag: splicing_instability
                        .genome_instability_splicing_flag[cell_id],
                }),
            },
            unspliced: unspliced.map(|u| JsonCellUnspliced {
                spliced_umis: u.spliced_umis[cell_id],
                unspliced_umis: u.unspliced_umis[cell_id],
                ambiguous_umis: u.ambiguous_umis[cell_id],
                unspliced_fraction: opt_f32(u.unspliced_fraction[cell_id]),
                unspliced_fraction_ci_low: opt_f32(u.unspliced_fraction_ci_low[cell_id]),
                unspliced_fraction_ci_high: opt_f32(u.unspliced_fraction_ci_high[cell_id]),
            }),
        });
    }

    let sis_stage = experimental.then(|| JsonSisStage {
        sis: opt_vec(&sis.sis),
        class: sis.class.iter().map(|c| class_str(*c)).collect(),
        p_missplicing: opt_vec(&sis.p_missplicing),
        p_imbalance: opt_vec(&sis.p_imbalance),
        p_entropy_z: opt_vec(&sis.p_entropy_z),
        p_entropy_abs: opt_vec(&sis.p_entropy_abs),
    });

    let isoform_stage = JsonIsoformStage {
        entropy: opt_vec(&isoform.entropy),
        dispersion: opt_vec(&isoform.dispersion),
        z_entropy: opt_vec(&isoform.z_entropy),
    };

    let missplicing_stage = JsonMissplicingStage {
        b_core: opt_vec(&missplicing.b_core),
        b_u12: opt_vec(&missplicing.b_u12),
        b_nmd: opt_vec(&missplicing.b_nmd),
        b_srhn: opt_vec(&missplicing.b_srhn),
        burden: opt_vec(&missplicing.burden),
        burden_star: opt_vec(&missplicing.burden_star),
    };

    let imbalance_stage = JsonImbalanceStage {
        z_u1: opt_vec(&imbalance.z_u1),
        z_u2: opt_vec(&imbalance.z_u2),
        z_sf3b: opt_vec(&imbalance.z_sf3b),
        z_srsf: opt_vec(&imbalance.z_srsf),
        z_hnrnp: opt_vec(&imbalance.z_hnrnp),
        z_u12: opt_vec(&imbalance.z_u12),
        z_nmd: opt_vec(&imbalance.z_nmd),
        axis_sr_hnrnp: opt_vec(&imbalance.axis_sr_hnrnp),
        axis_u2_u1: opt_vec(&imbalance.axis_u2_u1),
        axis_u12_major: opt_vec(&imbalance.axis_u12_major),
        axis_nmd: opt_vec(&imbalance.axis_nmd),
        imbalance: opt_vec(&imbalance.imbalance),
    };

    let splicing_noise_json = splicing_noise.map(|metrics| JsonSplicingNoise {
        noise_index: opt_f32(metrics.noise_index),
        per_geneset_noise: metrics
            .per_geneset_noise
            .iter()
            .map(|(geneset, noise)| JsonGenesetNoise {
                geneset: geneset.clone(),
                noise: opt_f32(*noise),
            })
            .collect(),
    });

    let cryptic_risk_json = cryptic_risk.map(|metrics| JsonCrypticRisk {
        cryptic_risk: opt_vec(&metrics.cryptic_risk),
        x_sr_hnrnp: opt_vec(&metrics.x_sr_hnrnp),
        x_entropy: opt_vec(&metrics.x_entropy),
        x_nmd: opt_vec(&metrics.x_nmd),
    });

    let collapse_json = collapse.map(|metrics| JsonCollapse {
        collapse_status: metrics
            .collapse_status
            .iter()
            .map(|s| collapse_str(*s))
            .collect(),
        core_suppression: metrics.core_suppression.clone(),
        high_imbalance: metrics.high_imbalance.clone(),
        low_sis: metrics.low_sis.clone(),
    });

    let timecourse_json = timecourse.map(|metrics| JsonTimecourse {
        trajectory: trajectory_str(metrics.trajectory),
        delta_sis: opt_vec(&metrics.delta_sis),
        delta_entropy: opt_vec(&metrics.delta_entropy),
        delta_imbalance: opt_vec(&metrics.delta_imbalance),
    });

    let splicing_instability_json =
        JsonSplicingInstabilityStage {
            panel_version: splicing_instability.panel_version,
            min_genes_per_panel_cell: splicing_instability.min_genes,
            conflict_panel_enabled: splicing_instability.conflict_panel_enabled,
            nmd_panel_enabled: splicing_instability.nmd_panel_enabled,
            splice_core: opt_vec(&splicing_instability.splice_core),
            rbp_core: opt_vec(&splicing_instability.rbp_core),
            rloop_resolve_core: opt_vec(&splicing_instability.rloop_resolve_core),
            conflict_risk_core: opt_vec(&splicing_instability.conflict_risk_core),
            nmd_core: opt_vec(&splicing_instability.nmd_core),
            composites: experimental.then(|| JsonSplicingInstabilityComposites {
                thresholds: JsonSplicingInstabilityThresholds {
                    splice_overload_high: SPLICE_OVERLOAD_HIGH_THRESHOLD,
                    rloop_risk_high: RLOOP_RISK_HIGH_THRESHOLD,
                    splicing_instability_high: SPLICING_INSTABILITY_HIGH_THRESHOLD,
                },
                sos: opt_vec(&splicing_instability.sos),
                rlr: opt_vec(&splicing_instability.rlr),
                sii: opt_vec(&splicing_instability.sii),
                splice_overload_high: splicing_instability.splice_overload_high.clone(),
                rloop_risk_high: splicing_instability.rloop_risk_high.clone(),
                splicing_instability_high: splicing_instability.splicing_instability_high.clone(),
                genome_instability_splicing_flag: splicing_instability
                    .genome_instability_splicing_flag
                    .clone(),
                global_stats: JsonSplicingInstabilityGlobalStats {
                    sos_p50: opt_f32(splicing_instability.global_stats.sos_p50),
                    sos_p90: opt_f32(splicing_instability.global_stats.sos_p90),
                    rlr_p50: opt_f32(splicing_instability.global_stats.rlr_p50),
                    rlr_p90: opt_f32(splicing_instability.global_stats.rlr_p90),
                    sii_p50: opt_f32(splicing_instability.global_stats.sii_p50),
                    sii_p90: opt_f32(splicing_instability.global_stats.sii_p90),
                },
            }),
            z_reference: JsonSplicingInstabilityZReference {
                splice_core: JsonRobustRef {
                    median: opt_f32(splicing_instability.z_reference.splice_core.median),
                    mad: opt_f32(splicing_instability.z_reference.splice_core.mad),
                },
                rbp_core: JsonRobustRef {
                    median: opt_f32(splicing_instability.z_reference.rbp_core.median),
                    mad: opt_f32(splicing_instability.z_reference.rbp_core.mad),
                },
                rloop_resolve_core: JsonRobustRef {
                    median: opt_f32(splicing_instability.z_reference.rloop_resolve_core.median),
                    mad: opt_f32(splicing_instability.z_reference.rloop_resolve_core.mad),
                },
                conflict_risk_core: splicing_instability
                    .z_reference
                    .conflict_risk_core
                    .as_ref()
                    .map(|r| JsonRobustRef {
                        median: opt_f32(r.median),
                        mad: opt_f32(r.mad),
                    }),
                nmd_core: splicing_instability.z_reference.nmd_core.as_ref().map(|r| {
                    JsonRobustRef {
                        median: opt_f32(r.median),
                        mad: opt_f32(r.mad),
                    }
                }),
            },
            missingness: JsonSplicingInstabilityMissingness {
                splice_core_nan_cells: splicing_instability.missingness.splice_core_nan_cells,
                rbp_core_nan_cells: splicing_instability.missingness.rbp_core_nan_cells,
                rloop_resolve_core_nan_cells: splicing_instability
                    .missingness
                    .rloop_resolve_core_nan_cells,
                conflict_risk_core_nan_cells: splicing_instability
                    .missingness
                    .conflict_risk_core_nan_cells,
                nmd_core_nan_cells: splicing_instability.missingness.nmd_core_nan_cells,
                sos_nan_cells: splicing_instability.missingness.sos_nan_cells,
                rlr_nan_cells: splicing_instability.missingness.rlr_nan_cells,
                sii_nan_cells: splicing_instability.missingness.sii_nan_cells,
                panel_coverage: splicing_instability
                    .missingness
                    .panel_coverage
                    .iter()
                    .map(|panel| JsonPanelCoverage {
                        panel_name: panel.panel_name,
                        genes_defined: panel.genes_defined,
                        genes_mapped: panel.genes_mapped,
                    })
                    .collect(),
            },
        };

    let output = JsonOutput {
        schema_version: JSON_SCHEMA_VERSION,
        tool: "kira-spliceqc",
        mode: "cell",
        n_cells,
        experimental_signatures: experimental,
        cells,
        sis: sis_stage,
        isoform: isoform_stage,
        missplicing: missplicing_stage,
        imbalance: imbalance_stage,
        splicing_noise: splicing_noise_json,
        cryptic_risk: cryptic_risk_json,
        collapse: collapse_json,
        timecourse: timecourse_json,
        splicing_instability: splicing_instability_json,
        cell_cycle_guardrail: None,
        input_levels: if unspliced.is_some() { vec!["L0", "L1"] } else { vec!["L0"] },
        unspliced: unspliced.map(|u| JsonUnsplicedStage {
            source: u.source.clone(),
            min_layer_umis: u.min_layer_umis,
            has_ambiguous: u.has_ambiguous,
            cells_without_layers: u.cells_without_layers,
            undefined_cells: u.undefined_cells,
            spliced_umis: u.spliced_umis.clone(),
            unspliced_umis: u.unspliced_umis.clone(),
            ambiguous_umis: u.ambiguous_umis.clone(),
            unspliced_fraction: opt_vec(&u.unspliced_fraction),
            unspliced_fraction_ci_low: opt_vec(&u.unspliced_fraction_ci_low),
            unspliced_fraction_ci_high: opt_vec(&u.unspliced_fraction_ci_high),
        }),
    };

    let file = std::fs::File::create(path).map_err(|e| InputError::io(path, e))?;
    let mut w = std::io::BufWriter::with_capacity(1 << 16, file);
    serde_json::to_writer(&mut w, &output)
        .map_err(|e| InputError::OutputSerialization(e.to_string()))?;
    use std::io::Write;
    w.flush().map_err(|e| InputError::io(path, e))?;
    Ok(())
}

fn opt_f32(value: f32) -> Option<f32> {
    if value.is_finite() { Some(value) } else { None }
}

fn opt_vec(values: &[f32]) -> Vec<Option<f32>> {
    values.iter().copied().map(opt_f32).collect()
}

fn class_str(class: SpliceIntegrityClass) -> &'static str {
    match class {
        SpliceIntegrityClass::Intact => "Intact",
        SpliceIntegrityClass::Stressed => "Stressed",
        SpliceIntegrityClass::Impaired => "Impaired",
        SpliceIntegrityClass::Broken => "Broken",
    }
}

fn collapse_str(status: SpliceosomeCollapseStatus) -> &'static str {
    match status {
        SpliceosomeCollapseStatus::NoCollapse => "NoCollapse",
        SpliceosomeCollapseStatus::Collapse => "Collapse",
        SpliceosomeCollapseStatus::Inconclusive => "Inconclusive",
    }
}

fn trajectory_str(class: SplicingTrajectoryClass) -> &'static str {
    match class {
        SplicingTrajectoryClass::Adaptive => "Adaptive",
        SplicingTrajectoryClass::Degenerative => "Degenerative",
        SplicingTrajectoryClass::Oscillatory => "Oscillatory",
        SplicingTrajectoryClass::Inconclusive => "Inconclusive",
    }
}
