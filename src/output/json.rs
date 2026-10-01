use std::path::Path;

use serde::Serialize;

use crate::input::error::InputError;
use crate::model::assembly_phase::AssemblyPhaseImbalanceMetrics;
use crate::model::cell_cycle::{CellCycleMetrics, CellCyclePhase};
use crate::model::cell_qc::CellQc;
use crate::model::collapse::{SpliceosomeCollapseMetrics, SpliceosomeCollapseStatus};
use crate::model::coupling::CouplingStressMetrics;
use crate::model::cryptic_risk::CrypticSplicingRiskMetrics;
use crate::model::exon_intron_bias::ExonIntronDefinitionMetrics;
use crate::model::imbalance::SpliceosomeImbalanceMetrics;
use crate::model::isoform_dispersion::IsoformDispersionMetrics;
use crate::model::junctions::JunctionMetrics;
use crate::model::missplicing::MissplicingMetrics;
use crate::model::sis::{SpliceIntegrityClass, SpliceIntegrityMetrics};
use crate::model::splicing_instability::{
    RLOOP_RISK_HIGH_THRESHOLD, SPLICE_OVERLOAD_HIGH_THRESHOLD, SPLICING_INSTABILITY_HIGH_THRESHOLD,
    SplicingInstabilityMetrics,
};
use crate::model::splicing_noise::SplicingNoiseMetrics;
use crate::model::timecourse::{SplicingTrajectoryClass, TimecourseSplicingMetrics};
use crate::model::intron_retention::IntronRetentionMetrics;
use crate::model::unspliced::UnsplicedMetrics;
use crate::output::provenance::{Provenance, input_levels};
use crate::reference::{MIN_STRATUM_CELLS, Strata};

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
    /// Cell-cycle confounder annotation (Tirosh 2016 scores, Seurat phase rule).
    cell_cycle: JsonCellCycleStage,
    /// Input levels available: `["L0"]` or `["L0", "L1"]`.
    input_levels: Vec<&'static str>,
    /// Tier A unspliced metrics (input level L1); absent without layers.
    #[serde(skip_serializing_if = "Option::is_none")]
    unspliced: Option<JsonUnsplicedStage>,
    /// Tier A intron retention (input level L1); absent without layers.
    #[serde(skip_serializing_if = "Option::is_none")]
    intron_retention: Option<JsonIntronRetentionStage>,
    /// Tier B junction metrics (input level L2); absent without a junction matrix.
    #[serde(skip_serializing_if = "Option::is_none")]
    junctions: Option<JsonJunctionStage>,
    /// Reference strata used for `_dev` metrics and outlier flags.
    reference: JsonReference,
    cell_qc: JsonCellQc,
    provenance: &'a Provenance,
}

#[derive(Serialize)]
struct JsonCellQc {
    min_counts: u64,
    min_genes: u64,
    doublet_column: Option<String>,
    n_low_depth: usize,
    n_doublet: usize,
    low_depth: Vec<bool>,
    doublet: Vec<bool>,
}

#[derive(Serialize)]
struct JsonCellJunctions {
    junction_umis: u64,
    unannotated_junction_fraction: Option<f32>,
    cryptic_3ss_umis: u64,
    cryptic_3ss_fraction: Option<f32>,
    cryptic_3ss_fraction_dev: Option<f32>,
    cryptic_3ss_high: bool,
    exon_skip_umis: u64,
    exon_skip_fraction: Option<f32>,
    exon_skip_fraction_dev: Option<f32>,
    exon_skip_high: bool,
    splice_site_shift: Option<f32>,
    splice_site_shift_dev: Option<f32>,
    splice_site_shift_high: bool,
    site_groups_used: u32,
}

#[derive(Serialize)]
struct JsonJunctionStage {
    source: String,
    n_junctions: usize,
    n_annotated: usize,
    n_cryptic_acceptor_junctions: usize,
    n_skip_junctions: usize,
    n_site_groups: usize,
    min_junction_umis: u64,
    min_ratio_umis: u64,
    cells_without_junctions: usize,
    undefined_cells: usize,
    junction_umis: Vec<u64>,
    unannotated_junction_fraction: Vec<Option<f32>>,
    cryptic_3ss_fraction: Vec<Option<f32>>,
    cryptic_3ss_fraction_dev: Vec<Option<f32>>,
    cryptic_3ss_high: Vec<bool>,
    exon_skip_fraction: Vec<Option<f32>>,
    exon_skip_fraction_dev: Vec<Option<f32>>,
    exon_skip_high: Vec<bool>,
    splice_site_shift: Vec<Option<f32>>,
    splice_site_shift_dev: Vec<Option<f32>>,
    splice_site_shift_high: Vec<bool>,
    cryptic_reference: Vec<JsonStratumStat>,
    skip_reference: Vec<JsonStratumStat>,
    shift_reference: Vec<JsonStratumStat>,
}

#[derive(Serialize)]
struct JsonCellIntronRetention {
    intron_retention_index: Option<f32>,
    intron_retention_index_dev: Option<f32>,
    ir_gene_dispersion: Option<f32>,
    ir_genes_used: u32,
    intron_retention_high: bool,
}

#[derive(Serialize)]
struct JsonIntronRetentionStage {
    min_gene_umis: u32,
    min_genes: usize,
    genes_with_reference: usize,
    undefined_cells: usize,
    intron_retention_index: Vec<Option<f32>>,
    intron_retention_index_dev: Vec<Option<f32>>,
    ir_gene_dispersion: Vec<Option<f32>>,
    ir_genes_used: Vec<u32>,
    intron_retention_high: Vec<bool>,
    reference: Vec<JsonStratumStat>,
}

#[derive(Serialize)]
struct JsonReference {
    mode: &'static str,
    column: Option<String>,
    n_strata: usize,
    min_stratum_cells: usize,
    folded_cells: usize,
    strata: Vec<JsonStratum>,
    /// Stratum label per cell (index into `strata`).
    labels: Vec<u32>,
}

#[derive(Serialize)]
struct JsonStratum {
    name: String,
    n_cells: usize,
}

#[derive(Serialize)]
struct JsonStratumStat {
    name: String,
    n_cells: usize,
    n_defined: usize,
    median: Option<f32>,
    mad: Option<f32>,
}

#[derive(Serialize)]
struct JsonCellUnspliced {
    spliced_umis: u64,
    unspliced_umis: u64,
    ambiguous_umis: u64,
    unspliced_fraction: Option<f32>,
    unspliced_fraction_ci_low: Option<f32>,
    unspliced_fraction_ci_high: Option<f32>,
    unspliced_fraction_dev: Option<f32>,
    nuclear_fraction_flag: bool,
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
    unspliced_fraction_dev: Vec<Option<f32>>,
    nuclear_fraction_flag: Vec<bool>,
    /// Per-stratum reference (median / MAD of `unspliced_fraction`).
    reference: Vec<JsonStratumStat>,
}

#[derive(Serialize)]
struct JsonCell<'a> {
    cell_id: usize,
    cell_name: &'a str,
    low_depth: bool,
    doublet: bool,
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
    #[serde(skip_serializing_if = "Option::is_none")]
    intron_retention: Option<JsonCellIntronRetention>,
    #[serde(skip_serializing_if = "Option::is_none")]
    junctions: Option<JsonCellJunctions>,
    cell_cycle: JsonCellCellCycle,
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
struct JsonCellCellCycle {
    s_score_expr: Option<f32>,
    g2m_score_expr: Option<f32>,
    phase: &'static str,
    cycling: bool,
}

#[derive(Serialize)]
struct JsonCellCycleStage {
    s_genes_mapped: usize,
    g2m_genes_mapped: usize,
    s_score_expr: Vec<Option<f32>>,
    g2m_score_expr: Vec<Option<f32>>,
    phase: Vec<&'static str>,
    cycling: Vec<bool>,
    phase_counts: JsonPhaseCounts,
}

#[derive(Serialize)]
struct JsonPhaseCounts {
    g1: usize,
    s: usize,
    g2m: usize,
    unknown: usize,
}

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
    intron_retention: Option<&IntronRetentionMetrics>,
    cell_cycle: &CellCycleMetrics,
    junctions: Option<&JunctionMetrics>,
    cell_qc: &CellQc,
    strata: &Strata,
    provenance: &Provenance,
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
            low_depth: cell_qc.low_depth[cell_id],
            doublet: cell_qc.doublet[cell_id],
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
                unspliced_fraction_dev: opt_f32(u.unspliced_fraction_dev[cell_id]),
                nuclear_fraction_flag: u.nuclear_fraction_flag[cell_id],
            }),
            intron_retention: intron_retention.map(|m| JsonCellIntronRetention {
                intron_retention_index: opt_f32(m.intron_retention_index[cell_id]),
                intron_retention_index_dev: opt_f32(m.intron_retention_index_dev[cell_id]),
                ir_gene_dispersion: opt_f32(m.ir_gene_dispersion[cell_id]),
                ir_genes_used: m.ir_genes_used[cell_id],
                intron_retention_high: m.intron_retention_high[cell_id],
            }),
            junctions: junctions.map(|m| JsonCellJunctions {
                junction_umis: m.junction_umis[cell_id],
                unannotated_junction_fraction: opt_f32(m.unannotated_junction_fraction[cell_id]),
                cryptic_3ss_umis: m.cryptic_3ss_umis[cell_id],
                cryptic_3ss_fraction: opt_f32(m.cryptic_3ss_fraction[cell_id]),
                cryptic_3ss_fraction_dev: opt_f32(m.cryptic_3ss_fraction_dev[cell_id]),
                cryptic_3ss_high: m.cryptic_3ss_high[cell_id],
                exon_skip_umis: m.exon_skip_umis[cell_id],
                exon_skip_fraction: opt_f32(m.exon_skip_fraction[cell_id]),
                exon_skip_fraction_dev: opt_f32(m.exon_skip_fraction_dev[cell_id]),
                exon_skip_high: m.exon_skip_high[cell_id],
                splice_site_shift: opt_f32(m.splice_site_shift[cell_id]),
                splice_site_shift_dev: opt_f32(m.splice_site_shift_dev[cell_id]),
                splice_site_shift_high: m.splice_site_shift_high[cell_id],
                site_groups_used: m.site_groups_used[cell_id],
            }),
            cell_cycle: JsonCellCellCycle {
                s_score_expr: opt_f32(cell_cycle.s_score[cell_id]),
                g2m_score_expr: opt_f32(cell_cycle.g2m_score[cell_id]),
                phase: cell_cycle.phase[cell_id].as_str(),
                cycling: cell_cycle.cycling[cell_id],
            },
        });
    }

    let phase_count = |p: CellCyclePhase| cell_cycle.phase.iter().filter(|q| **q == p).count();
    let cell_cycle_stage = JsonCellCycleStage {
        s_genes_mapped: cell_cycle.s_genes_mapped,
        g2m_genes_mapped: cell_cycle.g2m_genes_mapped,
        s_score_expr: opt_vec(&cell_cycle.s_score),
        g2m_score_expr: opt_vec(&cell_cycle.g2m_score),
        phase: cell_cycle.phase.iter().map(|p| p.as_str()).collect(),
        cycling: cell_cycle.cycling.clone(),
        phase_counts: JsonPhaseCounts {
            g1: phase_count(CellCyclePhase::G1),
            s: phase_count(CellCyclePhase::S),
            g2m: phase_count(CellCyclePhase::G2M),
            unknown: phase_count(CellCyclePhase::Unknown),
        },
    };

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
        cell_cycle: cell_cycle_stage,
        input_levels: input_levels(unspliced.is_some(), junctions.is_some()),
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
            unspliced_fraction_dev: opt_vec(&u.unspliced_fraction_dev),
            nuclear_fraction_flag: u.nuclear_fraction_flag.clone(),
            reference: u
                .reference
                .iter()
                .map(|s| JsonStratumStat {
                    name: s.name.clone(),
                    n_cells: s.n_cells,
                    n_defined: s.n_defined,
                    median: opt_f32(s.median),
                    mad: opt_f32(s.mad),
                })
                .collect(),
        }),
        intron_retention: intron_retention.map(|m| JsonIntronRetentionStage {
            min_gene_umis: m.min_gene_umis,
            min_genes: m.min_genes,
            genes_with_reference: m.genes_with_reference,
            undefined_cells: m.undefined_cells,
            intron_retention_index: opt_vec(&m.intron_retention_index),
            intron_retention_index_dev: opt_vec(&m.intron_retention_index_dev),
            ir_gene_dispersion: opt_vec(&m.ir_gene_dispersion),
            ir_genes_used: m.ir_genes_used.clone(),
            intron_retention_high: m.intron_retention_high.clone(),
            reference: m
                .reference
                .iter()
                .map(|s| JsonStratumStat {
                    name: s.name.clone(),
                    n_cells: s.n_cells,
                    n_defined: s.n_defined,
                    median: opt_f32(s.median),
                    mad: opt_f32(s.mad),
                })
                .collect(),
        }),
        junctions: junctions.map(|m| {
            let stats = |v: &[crate::reference::StratumStat]| -> Vec<JsonStratumStat> {
                v.iter()
                    .map(|s| JsonStratumStat {
                        name: s.name.clone(),
                        n_cells: s.n_cells,
                        n_defined: s.n_defined,
                        median: opt_f32(s.median),
                        mad: opt_f32(s.mad),
                    })
                    .collect()
            };
            JsonJunctionStage {
                source: m.source.clone(),
                n_junctions: m.n_junctions,
                n_annotated: m.n_annotated,
                n_cryptic_acceptor_junctions: m.n_cryptic_acceptor_junctions,
                n_skip_junctions: m.n_skip_junctions,
                n_site_groups: m.n_site_groups,
                min_junction_umis: m.min_junction_umis,
                min_ratio_umis: m.min_ratio_umis,
                cells_without_junctions: m.cells_without_junctions,
                undefined_cells: m.undefined_cells,
                junction_umis: m.junction_umis.clone(),
                unannotated_junction_fraction: opt_vec(&m.unannotated_junction_fraction),
                cryptic_3ss_fraction: opt_vec(&m.cryptic_3ss_fraction),
                cryptic_3ss_fraction_dev: opt_vec(&m.cryptic_3ss_fraction_dev),
                cryptic_3ss_high: m.cryptic_3ss_high.clone(),
                exon_skip_fraction: opt_vec(&m.exon_skip_fraction),
                exon_skip_fraction_dev: opt_vec(&m.exon_skip_fraction_dev),
                exon_skip_high: m.exon_skip_high.clone(),
                splice_site_shift: opt_vec(&m.splice_site_shift),
                splice_site_shift_dev: opt_vec(&m.splice_site_shift_dev),
                splice_site_shift_high: m.splice_site_shift_high.clone(),
                cryptic_reference: stats(&m.cryptic_reference),
                skip_reference: stats(&m.skip_reference),
                shift_reference: stats(&m.shift_reference),
            }
        }),
        reference: JsonReference {
            mode: strata.mode.as_str(),
            column: strata.column.clone(),
            n_strata: strata.n_strata(),
            min_stratum_cells: MIN_STRATUM_CELLS,
            folded_cells: strata.folded_cells,
            strata: strata
                .names
                .iter()
                .zip(strata.sizes())
                .map(|(name, n_cells)| JsonStratum {
                    name: name.clone(),
                    n_cells,
                })
                .collect(),
            labels: strata.labels.clone(),
        },
        cell_qc: JsonCellQc {
            min_counts: cell_qc.min_counts,
            min_genes: cell_qc.min_genes,
            doublet_column: cell_qc.doublet_column.clone(),
            n_low_depth: cell_qc.n_low_depth(),
            n_doublet: cell_qc.n_doublet(),
            low_depth: cell_qc.low_depth.clone(),
            doublet: cell_qc.doublet.clone(),
        },
        provenance,
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
