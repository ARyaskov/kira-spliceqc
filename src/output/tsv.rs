use std::path::Path;

use crate::input::error::InputError;
use crate::model::assembly_phase::AssemblyPhaseImbalanceMetrics;
use crate::model::coupling::CouplingStressMetrics;
use crate::model::exon_intron_bias::ExonIntronDefinitionMetrics;
use crate::model::imbalance::SpliceosomeImbalanceMetrics;
use crate::model::isoform_dispersion::IsoformDispersionMetrics;
use crate::model::missplicing::MissplicingMetrics;
use crate::model::sis::{SpliceIntegrityClass, SpliceIntegrityMetrics};
use crate::model::splicing_instability::SplicingInstabilityMetrics;

const HEADER: &str = "cell_id\tcell_name\tsis\tclass\tp_missplicing\tp_imbalance\tp_entropy_z\tp_entropy_abs\tiso_entropy\tiso_dispersion\tmissplicing_burden\timbalance\tcoupling_stress\texon_definition_bias\tea_imbalance\tb_imbalance\tcat_imbalance\tsplice_core\trbp_core\trloop_resolve_core\tconflict_risk_core\tnmd_core\tSOS\tRLR\tSII\tsplice_overload_high\trloop_risk_high\tsplicing_instability_high\tgenome_instability_splicing_flag";

pub fn write_tsv(
    path: &Path,
    cell_names: &[String],
    isoform: &IsoformDispersionMetrics,
    missplicing: &MissplicingMetrics,
    imbalance: &SpliceosomeImbalanceMetrics,
    sis: &SpliceIntegrityMetrics,
    splicing_instability: &SplicingInstabilityMetrics,
    coupling: &CouplingStressMetrics,
    exon_intron: &ExonIntronDefinitionMetrics,
    assembly: &AssemblyPhaseImbalanceMetrics,
) -> Result<(), InputError> {
    let n_cells = cell_names.len();
    let mut out = String::new();
    out.push_str(HEADER);
    out.push('\n');

    for cell_id in 0..n_cells {
        let row = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            cell_id,
            cell_names[cell_id],
            fmt_f32(sis.sis[cell_id]),
            class_str(sis.class[cell_id]),
            fmt_f32(sis.p_missplicing[cell_id]),
            fmt_f32(sis.p_imbalance[cell_id]),
            fmt_f32(sis.p_entropy_z[cell_id]),
            fmt_f32(sis.p_entropy_abs[cell_id]),
            fmt_f32(isoform.entropy[cell_id]),
            fmt_f32(isoform.dispersion[cell_id]),
            fmt_f32(missplicing.burden[cell_id]),
            fmt_f32(imbalance.imbalance[cell_id]),
            fmt_f32(coupling.coupling_stress[cell_id]),
            fmt_f32(exon_intron.exon_definition_bias[cell_id]),
            fmt_f32(assembly.ea_imbalance[cell_id]),
            fmt_f32(assembly.b_imbalance[cell_id]),
            fmt_f32(assembly.cat_imbalance[cell_id]),
            fmt_f32(splicing_instability.splice_core[cell_id]),
            fmt_f32(splicing_instability.rbp_core[cell_id]),
            fmt_f32(splicing_instability.rloop_resolve_core[cell_id]),
            fmt_f32(splicing_instability.conflict_risk_core[cell_id]),
            fmt_f32(splicing_instability.nmd_core[cell_id]),
            fmt_f32(splicing_instability.sos[cell_id]),
            fmt_f32(splicing_instability.rlr[cell_id]),
            fmt_f32(splicing_instability.sii[cell_id]),
            fmt_bool(splicing_instability.splice_overload_high[cell_id]),
            fmt_bool(splicing_instability.rloop_risk_high[cell_id]),
            fmt_bool(splicing_instability.splicing_instability_high[cell_id]),
            fmt_bool(splicing_instability.genome_instability_splicing_flag[cell_id]),
        );
        out.push_str(&row);
    }

    std::fs::write(path, out).map_err(|e| InputError::io(path, e))?;
    Ok(())
}

pub fn header() -> &'static str {
    HEADER
}

fn fmt_f32(value: f32) -> String {
    if value.is_finite() {
        value.to_string()
    } else {
        String::new()
    }
}

fn class_str(class: SpliceIntegrityClass) -> &'static str {
    match class {
        SpliceIntegrityClass::Intact => "Intact",
        SpliceIntegrityClass::Stressed => "Stressed",
        SpliceIntegrityClass::Impaired => "Impaired",
        SpliceIntegrityClass::Broken => "Broken",
    }
}

fn fmt_bool(value: bool) -> &'static str {
    if value { "true" } else { "false" }
}
