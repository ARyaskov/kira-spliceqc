use std::path::Path;
use std::time::Instant;

use tracing::{debug, info};

use crate::input::error::InputError;
use crate::model::assembly_phase::AssemblyPhaseImbalanceMetrics;
use crate::model::collapse::SpliceosomeCollapseMetrics;
use crate::model::coupling::CouplingStressMetrics;
use crate::model::cryptic_risk::CrypticSplicingRiskMetrics;
use crate::model::exon_intron_bias::ExonIntronDefinitionMetrics;
use crate::model::imbalance::SpliceosomeImbalanceMetrics;
use crate::model::isoform_dispersion::IsoformDispersionMetrics;
use crate::model::missplicing::MissplicingMetrics;
use crate::model::sis::SpliceIntegrityMetrics;
use crate::model::splicing_noise::SplicingNoiseMetrics;
use crate::model::timecourse::TimecourseSplicingMetrics;
use crate::output::{json, summary, tsv};

pub struct OutputOptions {
    pub json: bool,
    pub tsv: bool,
}

pub fn run_stage7(
    out_dir: &Path,
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
    options: OutputOptions,
) -> Result<String, InputError> {
    validate_lengths(
        cell_names,
        isoform,
        missplicing,
        imbalance,
        sis,
        coupling,
        exon_intron,
        assembly,
        cryptic_risk,
        collapse,
    )?;
    std::fs::create_dir_all(out_dir).map_err(|e| InputError::io(out_dir, e))?;

    let write_json = options.json || (!options.json && !options.tsv);
    let write_tsv = options.tsv || (!options.json && !options.tsv);

    let start = Instant::now();

    if write_json {
        let path = out_dir.join("spliceqc.json");
        json::write_json(
            &path,
            cell_names,
            isoform,
            missplicing,
            imbalance,
            sis,
            coupling,
            exon_intron,
            assembly,
            splicing_noise,
            cryptic_risk,
            collapse,
            timecourse,
        )?;
        info!("wrote {}", path.display());
    }

    let coupling_default;
    let exon_intron_default;
    let assembly_default;

    let coupling_ref = match coupling {
        Some(value) => value,
        None => {
            coupling_default = CouplingStressMetrics {
                coupling_stress: vec![f32::NAN; cell_names.len()],
            };
            &coupling_default
        }
    };

    let exon_intron_ref = match exon_intron {
        Some(value) => value,
        None => {
            exon_intron_default = ExonIntronDefinitionMetrics {
                exon_definition_bias: vec![f32::NAN; cell_names.len()],
                z_srsf: vec![f32::NAN; cell_names.len()],
                z_u2af: vec![f32::NAN; cell_names.len()],
                z_hnrnp: vec![f32::NAN; cell_names.len()],
            };
            &exon_intron_default
        }
    };

    let assembly_ref = match assembly {
        Some(value) => value,
        None => {
            assembly_default = AssemblyPhaseImbalanceMetrics {
                z_ea: vec![f32::NAN; cell_names.len()],
                z_b: vec![f32::NAN; cell_names.len()],
                z_cat: vec![f32::NAN; cell_names.len()],
                ea_imbalance: vec![f32::NAN; cell_names.len()],
                b_imbalance: vec![f32::NAN; cell_names.len()],
                cat_imbalance: vec![f32::NAN; cell_names.len()],
            };
            &assembly_default
        }
    };

    if write_tsv {
        let path = out_dir.join("spliceqc.tsv");
        tsv::write_tsv(
            &path,
            cell_names,
            isoform,
            missplicing,
            imbalance,
            sis,
            coupling_ref,
            exon_intron_ref,
            assembly_ref,
        )?;
        info!("wrote {}", path.display());
    }

    let summary_text = summary::format_summary(sis, cryptic_risk, collapse, None);
    debug!(
        elapsed_ms = start.elapsed().as_millis(),
        "output serialization complete"
    );

    Ok(summary_text)
}

fn validate_lengths(
    cell_names: &[String],
    isoform: &IsoformDispersionMetrics,
    missplicing: &MissplicingMetrics,
    imbalance: &SpliceosomeImbalanceMetrics,
    sis: &SpliceIntegrityMetrics,
    coupling: Option<&CouplingStressMetrics>,
    exon_intron: Option<&ExonIntronDefinitionMetrics>,
    assembly: Option<&AssemblyPhaseImbalanceMetrics>,
    cryptic_risk: Option<&CrypticSplicingRiskMetrics>,
    collapse: Option<&SpliceosomeCollapseMetrics>,
) -> Result<(), InputError> {
    let n = cell_names.len();
    if isoform.entropy.len() != n || isoform.z_entropy.len() != n || isoform.dispersion.len() != n {
        return Err(InputError::LengthMismatch(
            "isoform length mismatch".to_string(),
        ));
    }
    if missplicing.burden.len() != n || missplicing.burden_star.len() != n {
        return Err(InputError::LengthMismatch(
            "missplicing length mismatch".to_string(),
        ));
    }
    if imbalance.imbalance.len() != n {
        return Err(InputError::LengthMismatch(
            "imbalance length mismatch".to_string(),
        ));
    }
    if sis.sis.len() != n || sis.class.len() != n {
        return Err(InputError::LengthMismatch(
            "sis length mismatch".to_string(),
        ));
    }
    if let Some(coupling) = coupling {
        if coupling.coupling_stress.len() != n {
            return Err(InputError::LengthMismatch(
                "coupling length mismatch".to_string(),
            ));
        }
    }
    if let Some(exon_intron) = exon_intron {
        if exon_intron.exon_definition_bias.len() != n {
            return Err(InputError::LengthMismatch(
                "exon/intron length mismatch".to_string(),
            ));
        }
    }
    if let Some(assembly) = assembly {
        if assembly.ea_imbalance.len() != n {
            return Err(InputError::LengthMismatch(
                "assembly phase length mismatch".to_string(),
            ));
        }
    }
    if let Some(cryptic) = cryptic_risk {
        if cryptic.cryptic_risk.len() != n
            || cryptic.x_sr_hnrnp.len() != n
            || cryptic.x_entropy.len() != n
            || cryptic.x_nmd.len() != n
        {
            return Err(InputError::LengthMismatch(
                "cryptic risk length mismatch".to_string(),
            ));
        }
    }
    if let Some(collapse) = collapse {
        if collapse.collapse_status.len() != n
            || collapse.core_suppression.len() != n
            || collapse.high_imbalance.len() != n
            || collapse.low_sis.len() != n
        {
            return Err(InputError::LengthMismatch(
                "collapse length mismatch".to_string(),
            ));
        }
    }
    Ok(())
}
