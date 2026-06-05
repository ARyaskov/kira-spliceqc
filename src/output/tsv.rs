use std::fs::File;
use std::io::{BufWriter, Write};
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
    let file = File::create(path).map_err(|e| InputError::io(path, e))?;
    let mut w = BufWriter::with_capacity(1 << 16, file);

    writeln!(w, "{}", HEADER).map_err(|e| InputError::io(path, e))?;

    let n_cells = cell_names.len();
    let mut buf = ryu_buf();
    for cell_id in 0..n_cells {
        write!(w, "{}\t{}\t", cell_id, cell_names[cell_id])
            .map_err(|e| InputError::io(path, e))?;
        write_f32(&mut w, sis.sis[cell_id], &mut buf, path)?;
        w.write_all(b"\t").map_err(|e| InputError::io(path, e))?;
        write!(w, "{}", class_str(sis.class[cell_id]))
            .map_err(|e| InputError::io(path, e))?;
        for v in [
            sis.p_missplicing[cell_id],
            sis.p_imbalance[cell_id],
            sis.p_entropy_z[cell_id],
            sis.p_entropy_abs[cell_id],
            isoform.entropy[cell_id],
            isoform.dispersion[cell_id],
            missplicing.burden[cell_id],
            imbalance.imbalance[cell_id],
            coupling.coupling_stress[cell_id],
            exon_intron.exon_definition_bias[cell_id],
            assembly.ea_imbalance[cell_id],
            assembly.b_imbalance[cell_id],
            assembly.cat_imbalance[cell_id],
            splicing_instability.splice_core[cell_id],
            splicing_instability.rbp_core[cell_id],
            splicing_instability.rloop_resolve_core[cell_id],
            splicing_instability.conflict_risk_core[cell_id],
            splicing_instability.nmd_core[cell_id],
            splicing_instability.sos[cell_id],
            splicing_instability.rlr[cell_id],
            splicing_instability.sii[cell_id],
        ] {
            w.write_all(b"\t").map_err(|e| InputError::io(path, e))?;
            write_f32(&mut w, v, &mut buf, path)?;
        }
        for v in [
            splicing_instability.splice_overload_high[cell_id],
            splicing_instability.rloop_risk_high[cell_id],
            splicing_instability.splicing_instability_high[cell_id],
            splicing_instability.genome_instability_splicing_flag[cell_id],
        ] {
            w.write_all(b"\t").map_err(|e| InputError::io(path, e))?;
            w.write_all(if v { b"true" } else { b"false" })
                .map_err(|e| InputError::io(path, e))?;
        }
        w.write_all(b"\n").map_err(|e| InputError::io(path, e))?;
    }

    w.flush().map_err(|e| InputError::io(path, e))?;
    Ok(())
}

pub fn header() -> &'static str {
    HEADER
}

#[inline]
fn class_str(class: SpliceIntegrityClass) -> &'static str {
    match class {
        SpliceIntegrityClass::Intact => "Intact",
        SpliceIntegrityClass::Stressed => "Stressed",
        SpliceIntegrityClass::Impaired => "Impaired",
        SpliceIntegrityClass::Broken => "Broken",
    }
}

/// Scratch buffer for float→ascii formatting. Reused across the row loop to
/// avoid per-field String allocations (was ~29 × n_cells `format!()` calls).
#[inline]
fn ryu_buf() -> String {
    String::with_capacity(32)
}

#[inline]
fn write_f32<W: Write>(
    w: &mut W,
    value: f32,
    buf: &mut String,
    path: &Path,
) -> Result<(), InputError> {
    use std::fmt::Write as _;
    buf.clear();
    if value.is_finite() {
        // Preserves f32 round-trip; matches `value.to_string()` from previous impl
        // to keep the deterministic_output_hash test stable.
        let _ = write!(buf, "{}", value);
    }
    w.write_all(buf.as_bytes()).map_err(|e| InputError::io(path, e))?;
    Ok(())
}
