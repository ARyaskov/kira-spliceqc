use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::input::error::InputError;
use crate::model::assembly_phase::AssemblyPhaseImbalanceMetrics;
use crate::model::cell_cycle::CellCycleMetrics;
use crate::model::cell_qc::CellQc;
use crate::model::coupling::CouplingStressMetrics;
use crate::model::exon_intron_bias::ExonIntronDefinitionMetrics;
use crate::model::imbalance::SpliceosomeImbalanceMetrics;
use crate::model::intron_retention::IntronRetentionMetrics;
use crate::model::isoform_dispersion::IsoformDispersionMetrics;
use crate::model::junctions::JunctionMetrics;
use crate::model::missplicing::MissplicingMetrics;
use crate::model::sis::{SpliceIntegrityClass, SpliceIntegrityMetrics};
use crate::model::splicing_instability::SplicingInstabilityMetrics;
use crate::model::unspliced::UnsplicedMetrics;

/// Column naming: every metric derived purely from panel expression carries
/// the `_expr` suffix (it is an expression signature, not a measurement of
/// splicing). Composite indices (`sis`, `class`, `p_*`, `SOS`, `RLR`, `SII`)
/// and their flags are experimental (see METRICS.md) and are only written
/// when `experimental` is set.
///
/// `(name, experimental)` in output order. `cell_values` must produce values
/// in exactly this order.
const COLUMNS: &[(&str, bool)] = &[
    ("cell_id", false),
    ("cell_name", false),
    ("low_depth", false),
    ("doublet", false),
    ("sis", true),
    ("class", true),
    ("p_missplicing", true),
    ("p_imbalance", true),
    ("p_entropy_z", true),
    ("p_entropy_abs", true),
    ("regulator_entropy_expr", false),
    ("regulator_dispersion_expr", false),
    ("missplicing_burden_expr", false),
    ("spliceosome_imbalance_expr", false),
    ("coupling_stress_expr", false),
    ("exon_definition_bias_expr", false),
    ("ea_phase_imbalance_expr", false),
    ("b_phase_imbalance_expr", false),
    ("catalytic_phase_imbalance_expr", false),
    ("spliceosome_core_expr", false),
    ("splicing_rbp_expr", false),
    ("rloop_resolution_expr", false),
    ("conflict_risk_expr", false),
    ("nmd_factor_expr", false),
    // Tier A (input level L1); empty when no layers were loaded.
    ("spliced_umis", false),
    ("unspliced_umis", false),
    ("ambiguous_umis", false),
    ("unspliced_fraction", false),
    ("unspliced_fraction_ci_low", false),
    ("unspliced_fraction_ci_high", false),
    ("unspliced_fraction_dev", false),
    ("nuclear_fraction_flag", false),
    ("intron_retention_index", false),
    ("intron_retention_index_dev", false),
    ("ir_gene_dispersion", false),
    ("ir_genes_used", false),
    ("intron_retention_high", false),
    // Tier B (input level L2); empty without a junction matrix.
    ("junction_umis", false),
    ("unannotated_junction_fraction", false),
    ("cryptic_3ss_umis", false),
    ("cryptic_3ss_fraction", false),
    ("cryptic_3ss_fraction_dev", false),
    ("cryptic_3ss_high", false),
    ("exon_skip_umis", false),
    ("exon_skip_fraction", false),
    ("exon_skip_fraction_dev", false),
    ("exon_skip_high", false),
    ("splice_site_shift", false),
    ("splice_site_shift_dev", false),
    ("splice_site_shift_high", false),
    ("site_groups_used", false),
    // Cell-cycle confounder annotation (Tirosh 2016 / Seurat rule).
    ("s_score_expr", false),
    ("g2m_score_expr", false),
    ("cell_cycle_phase", false),
    ("cycling", false),
    ("SOS", true),
    ("RLR", true),
    ("SII", true),
    ("splice_overload_high", true),
    ("rloop_risk_high", true),
    ("splicing_instability_high", true),
    ("genome_instability_splicing_flag", true),
];

enum Value<'a> {
    Index(usize),
    Str(&'a str),
    F32(f32),
    Bool(bool),
    /// Optional integer count (empty field when None).
    OptU64(Option<u64>),
    /// Optional flag (empty field when None).
    OptBool(Option<bool>),
}

#[allow(clippy::too_many_arguments)]
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
    unspliced: Option<&UnsplicedMetrics>,
    intron_retention: Option<&IntronRetentionMetrics>,
    cell_cycle: &CellCycleMetrics,
    junctions: Option<&JunctionMetrics>,
    cell_qc: &CellQc,
    experimental: bool,
) -> Result<(), InputError> {
    let file = File::create(path).map_err(|e| InputError::io(path, e))?;
    let mut w = BufWriter::with_capacity(1 << 16, file);

    writeln!(w, "{}", header(experimental)).map_err(|e| InputError::io(path, e))?;

    let n_cells = cell_names.len();
    let mut buf = String::with_capacity(32);
    let mut values: Vec<Value> = Vec::with_capacity(COLUMNS.len());
    for cell_id in 0..n_cells {
        values.clear();
        cell_values(
            &mut values,
            cell_id,
            cell_names,
            isoform,
            missplicing,
            imbalance,
            sis,
            splicing_instability,
            coupling,
            exon_intron,
            assembly,
            unspliced,
            intron_retention,
            cell_cycle,
            junctions,
            cell_qc,
        );
        debug_assert_eq!(values.len(), COLUMNS.len());

        let mut first = true;
        for ((_, is_experimental), value) in COLUMNS.iter().zip(values.iter()) {
            if *is_experimental && !experimental {
                continue;
            }
            if !first {
                w.write_all(b"\t").map_err(|e| InputError::io(path, e))?;
            }
            first = false;
            write_value(&mut w, value, &mut buf, path)?;
        }
        w.write_all(b"\n").map_err(|e| InputError::io(path, e))?;
    }

    w.flush().map_err(|e| InputError::io(path, e))?;
    Ok(())
}

/// Header line for the given mode (composite columns only when `experimental`).
pub fn header(experimental: bool) -> String {
    COLUMNS
        .iter()
        .filter(|(_, is_experimental)| experimental || !*is_experimental)
        .map(|(name, _)| *name)
        .collect::<Vec<_>>()
        .join("\t")
}

#[allow(clippy::too_many_arguments)]
fn cell_values<'a>(
    out: &mut Vec<Value<'a>>,
    cell_id: usize,
    cell_names: &'a [String],
    isoform: &IsoformDispersionMetrics,
    missplicing: &MissplicingMetrics,
    imbalance: &SpliceosomeImbalanceMetrics,
    sis: &SpliceIntegrityMetrics,
    si: &SplicingInstabilityMetrics,
    coupling: &CouplingStressMetrics,
    exon_intron: &ExonIntronDefinitionMetrics,
    assembly: &AssemblyPhaseImbalanceMetrics,
    unspliced: Option<&UnsplicedMetrics>,
    ir: Option<&IntronRetentionMetrics>,
    cc: &CellCycleMetrics,
    jn: Option<&JunctionMetrics>,
    qc: &CellQc,
) {
    out.push(Value::Index(cell_id));
    out.push(Value::Str(&cell_names[cell_id]));
    out.push(Value::Bool(qc.low_depth[cell_id]));
    out.push(Value::Bool(qc.doublet[cell_id]));
    out.push(Value::F32(sis.sis[cell_id]));
    out.push(Value::Str(class_str(sis.class[cell_id])));
    out.push(Value::F32(sis.p_missplicing[cell_id]));
    out.push(Value::F32(sis.p_imbalance[cell_id]));
    out.push(Value::F32(sis.p_entropy_z[cell_id]));
    out.push(Value::F32(sis.p_entropy_abs[cell_id]));
    out.push(Value::F32(isoform.entropy[cell_id]));
    out.push(Value::F32(isoform.dispersion[cell_id]));
    out.push(Value::F32(missplicing.burden[cell_id]));
    out.push(Value::F32(imbalance.imbalance[cell_id]));
    out.push(Value::F32(coupling.coupling_stress[cell_id]));
    out.push(Value::F32(exon_intron.exon_definition_bias[cell_id]));
    out.push(Value::F32(assembly.ea_imbalance[cell_id]));
    out.push(Value::F32(assembly.b_imbalance[cell_id]));
    out.push(Value::F32(assembly.cat_imbalance[cell_id]));
    out.push(Value::F32(si.splice_core[cell_id]));
    out.push(Value::F32(si.rbp_core[cell_id]));
    out.push(Value::F32(si.rloop_resolve_core[cell_id]));
    out.push(Value::F32(si.conflict_risk_core[cell_id]));
    out.push(Value::F32(si.nmd_core[cell_id]));
    out.push(Value::OptU64(unspliced.map(|u| u.spliced_umis[cell_id])));
    out.push(Value::OptU64(unspliced.map(|u| u.unspliced_umis[cell_id])));
    out.push(Value::OptU64(unspliced.map(|u| u.ambiguous_umis[cell_id])));
    out.push(Value::F32(
        unspliced.map_or(f32::NAN, |u| u.unspliced_fraction[cell_id]),
    ));
    out.push(Value::F32(
        unspliced.map_or(f32::NAN, |u| u.unspliced_fraction_ci_low[cell_id]),
    ));
    out.push(Value::F32(
        unspliced.map_or(f32::NAN, |u| u.unspliced_fraction_ci_high[cell_id]),
    ));
    out.push(Value::F32(
        unspliced.map_or(f32::NAN, |u| u.unspliced_fraction_dev[cell_id]),
    ));
    out.push(Value::OptBool(
        unspliced.map(|u| u.nuclear_fraction_flag[cell_id]),
    ));
    out.push(Value::F32(
        ir.map_or(f32::NAN, |m| m.intron_retention_index[cell_id]),
    ));
    out.push(Value::F32(
        ir.map_or(f32::NAN, |m| m.intron_retention_index_dev[cell_id]),
    ));
    out.push(Value::F32(
        ir.map_or(f32::NAN, |m| m.ir_gene_dispersion[cell_id]),
    ));
    out.push(Value::OptU64(ir.map(|m| m.ir_genes_used[cell_id] as u64)));
    out.push(Value::OptBool(ir.map(|m| m.intron_retention_high[cell_id])));
    out.push(Value::OptU64(jn.map(|m| m.junction_umis[cell_id])));
    out.push(Value::F32(
        jn.map_or(f32::NAN, |m| m.unannotated_junction_fraction[cell_id]),
    ));
    out.push(Value::OptU64(jn.map(|m| m.cryptic_3ss_umis[cell_id])));
    out.push(Value::F32(
        jn.map_or(f32::NAN, |m| m.cryptic_3ss_fraction[cell_id]),
    ));
    out.push(Value::F32(
        jn.map_or(f32::NAN, |m| m.cryptic_3ss_fraction_dev[cell_id]),
    ));
    out.push(Value::OptBool(jn.map(|m| m.cryptic_3ss_high[cell_id])));
    out.push(Value::OptU64(jn.map(|m| m.exon_skip_umis[cell_id])));
    out.push(Value::F32(
        jn.map_or(f32::NAN, |m| m.exon_skip_fraction[cell_id]),
    ));
    out.push(Value::F32(
        jn.map_or(f32::NAN, |m| m.exon_skip_fraction_dev[cell_id]),
    ));
    out.push(Value::OptBool(jn.map(|m| m.exon_skip_high[cell_id])));
    out.push(Value::F32(
        jn.map_or(f32::NAN, |m| m.splice_site_shift[cell_id]),
    ));
    out.push(Value::F32(
        jn.map_or(f32::NAN, |m| m.splice_site_shift_dev[cell_id]),
    ));
    out.push(Value::OptBool(
        jn.map(|m| m.splice_site_shift_high[cell_id]),
    ));
    out.push(Value::OptU64(
        jn.map(|m| m.site_groups_used[cell_id] as u64),
    ));
    out.push(Value::F32(cc.s_score[cell_id]));
    out.push(Value::F32(cc.g2m_score[cell_id]));
    out.push(Value::Str(cc.phase[cell_id].as_str()));
    out.push(Value::Bool(cc.cycling[cell_id]));
    out.push(Value::F32(si.sos[cell_id]));
    out.push(Value::F32(si.rlr[cell_id]));
    out.push(Value::F32(si.sii[cell_id]));
    out.push(Value::Bool(si.splice_overload_high[cell_id]));
    out.push(Value::Bool(si.rloop_risk_high[cell_id]));
    out.push(Value::Bool(si.splicing_instability_high[cell_id]));
    out.push(Value::Bool(si.genome_instability_splicing_flag[cell_id]));
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

#[inline]
fn write_value<W: Write>(
    w: &mut W,
    value: &Value<'_>,
    buf: &mut String,
    path: &Path,
) -> Result<(), InputError> {
    use std::fmt::Write as _;
    buf.clear();
    match value {
        Value::Index(i) => {
            let _ = write!(buf, "{i}");
        }
        Value::Str(s) => buf.push_str(s),
        // Non-finite values are written as empty fields. `{}` on f32 is the
        // shortest round-trip representation, kept for output determinism.
        Value::F32(v) => {
            if v.is_finite() {
                let _ = write!(buf, "{v}");
            }
        }
        Value::Bool(b) => buf.push_str(if *b { "true" } else { "false" }),
        Value::OptU64(Some(v)) => {
            let _ = write!(buf, "{v}");
        }
        Value::OptU64(None) => {}
        Value::OptBool(Some(b)) => buf.push_str(if *b { "true" } else { "false" }),
        Value::OptBool(None) => {}
    }
    w.write_all(buf.as_bytes())
        .map_err(|e| InputError::io(path, e))
}
