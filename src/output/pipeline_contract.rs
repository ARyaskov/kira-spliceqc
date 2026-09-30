use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use rayon::prelude::*;
use serde::Serialize;
use tracing::info;

use crate::expression::ExpressionMatrix;
use crate::genesets::{Geneset, GenesetCatalog};
use crate::input::InputDescriptor;
use crate::input::error::InputError;
use crate::model::coupling::CouplingStressMetrics;
use crate::model::imbalance::SpliceosomeImbalanceMetrics;
use crate::model::missplicing::MissplicingMetrics;
use crate::model::sis::SpliceIntegrityMetrics;
use crate::model::intron_retention::IntronRetentionMetrics;
use crate::model::unspliced::UnsplicedMetrics;
use crate::reference::{MIN_STRATUM_CELLS, Strata};
use crate::model::splicing_instability::{
    RLOOP_RISK_HIGH_THRESHOLD, SPLICE_OVERLOAD_HIGH_THRESHOLD, SPLICING_INSTABILITY_HIGH_THRESHOLD,
    SplicingInstabilityMetrics,
};
use crate::stats::robust::quantile_f64;

const PIPELINE_DIR: &str = "kira-spliceqc";

/// Version of the spliceqc.tsv / summary.json contract consumed by kira-organelle.
pub const PIPELINE_CONTRACT_VERSION: &str = "0.3";

const SPLICEQC_HEADER: &str = "barcode\tsample\tcondition\tspecies\tlibsize\tnnz\texpressed_genes\tsplice_fidelity_index\tintron_retention_rate\texon_skipping_rate\talt_splice_burden\tsplice_junction_noise\tstress_splicing_index\tregime\tflags\tconfidence";

const REGIMES: [&str; 6] = [
    "HighFidelitySplicing",
    "RegulatedAlternativeSplicing",
    "StressInducedSplicing",
    "SpliceNoiseDominant",
    "SplicingCollapse",
    "Unclassified",
];

#[derive(Debug, Clone)]
pub struct PipelineCellRow {
    pub barcode: String,
    pub sample: String,
    pub condition: String,
    pub species: String,
    pub libsize: u64,
    pub nnz: u64,
    pub expressed_genes: u64,
    pub splice_fidelity_index: f64,
    pub intron_retention_rate: f64,
    pub exon_skipping_rate: f64,
    pub alt_splice_burden: f64,
    pub splice_junction_noise: f64,
    pub stress_splicing_index: f64,
    pub regime: &'static str,
    pub flags: String,
    pub confidence: f64,
}

#[derive(Serialize)]
struct SummaryJson {
    tool: ToolJson,
    input: InputJson,
    distributions: DistributionsJson,
    regimes: RegimesJson,
    qc: QcJson,
    splicing_instability: SplicingInstabilitySummaryJson,
    /// Tier A summary (input level L1); `null` without layers.
    unspliced: Option<UnsplicedSummaryJson>,
    /// Tier A intron retention summary; `null` without layers.
    intron_retention: Option<IntronRetentionSummaryJson>,
    /// Reference strata used for deviations and flags.
    reference: ReferenceJson,
}

#[derive(Serialize)]
struct IntronRetentionSummaryJson {
    min_gene_umis: u32,
    min_genes: usize,
    genes_with_reference: usize,
    n_defined_cells: usize,
    median: Option<f64>,
    p10: Option<f64>,
    p90: Option<f64>,
    high_flag_fraction: f64,
    strata: Vec<StratumStatJson>,
}

#[derive(Serialize)]
struct ToolJson {
    name: &'static str,
    version: &'static str,
    simd: &'static str,
    /// The contract metrics are built on experimental composite signatures.
    experimental_signatures: bool,
}

#[derive(Serialize)]
struct InputJson {
    n_cells: usize,
    species: &'static str,
    /// Input levels available: `["L0"]` or `["L0", "L1"]`.
    levels: Vec<&'static str>,
}

#[derive(Serialize)]
struct UnsplicedSummaryJson {
    source: String,
    min_layer_umis: u64,
    n_defined_cells: usize,
    cells_without_layers: usize,
    median: Option<f64>,
    p10: Option<f64>,
    p90: Option<f64>,
    nuclear_fraction_flag_fraction: f64,
    /// Per-stratum reference of the unspliced fraction.
    strata: Vec<StratumStatJson>,
}

#[derive(Serialize)]
struct StratumStatJson {
    name: String,
    n_cells: usize,
    n_defined: usize,
    median: Option<f32>,
    mad: Option<f32>,
}

#[derive(Serialize)]
struct ReferenceJson {
    mode: &'static str,
    column: Option<String>,
    n_strata: usize,
    min_stratum_cells: usize,
    folded_cells: usize,
    strata: Vec<StratumSizeJson>,
}

#[derive(Serialize)]
struct StratumSizeJson {
    name: String,
    n_cells: usize,
}

#[derive(Serialize)]
struct DistributionJson {
    median: f64,
    p90: f64,
    p99: f64,
}

#[derive(Serialize)]
struct DistributionsJson {
    splice_fidelity_index: DistributionJson,
    stress_splicing_index: DistributionJson,
}

#[derive(Serialize)]
struct RegimesJson {
    counts: BTreeMap<&'static str, u64>,
    fractions: BTreeMap<&'static str, f64>,
}

#[derive(Serialize)]
struct QcJson {
    low_confidence_fraction: f64,
    high_splice_noise_fraction: f64,
}

#[derive(Serialize)]
struct SplicingInstabilityThresholdsJson {
    splice_overload_high: f32,
    rloop_risk_high: f32,
    splicing_instability_high: f32,
}

#[derive(Serialize)]
struct SplicingInstabilityRobustRefJson {
    median: Option<f32>,
    mad: Option<f32>,
}

#[derive(Serialize)]
struct SplicingInstabilityZRefJson {
    splice_core: SplicingInstabilityRobustRefJson,
    rbp_core: SplicingInstabilityRobustRefJson,
    rloop_resolve_core: SplicingInstabilityRobustRefJson,
    conflict_risk_core: Option<SplicingInstabilityRobustRefJson>,
    nmd_core: Option<SplicingInstabilityRobustRefJson>,
}

#[derive(Serialize)]
struct SplicingInstabilityGlobalStatsJson {
    sos_p50: Option<f32>,
    sos_p90: Option<f32>,
    rlr_p50: Option<f32>,
    rlr_p90: Option<f32>,
    sii_p50: Option<f32>,
    sii_p90: Option<f32>,
}

#[derive(Serialize)]
struct SplicingInstabilityClusterStatsJson {
    cluster_id: String,
    n_cells: usize,
    sos_median: Option<f32>,
    sos_p10: Option<f32>,
    sos_p90: Option<f32>,
    rlr_median: Option<f32>,
    rlr_p10: Option<f32>,
    rlr_p90: Option<f32>,
    sii_median: Option<f32>,
    sii_p10: Option<f32>,
    sii_p90: Option<f32>,
    splice_overload_high_fraction: f32,
    rloop_risk_high_fraction: f32,
    splicing_instability_high_fraction: f32,
    genome_instability_splicing_flag_fraction: f32,
}

#[derive(Serialize)]
struct SplicingInstabilityTopClusterJson {
    by_sii_median: Vec<String>,
    by_genome_instability_splicing_flag_fraction: Vec<String>,
}

#[derive(Serialize)]
struct SplicingInstabilityPanelCoverageJson {
    panel_name: &'static str,
    genes_defined: usize,
    genes_mapped: usize,
}

#[derive(Serialize)]
struct SplicingInstabilityMissingnessJson {
    splice_core_nan_cells: usize,
    rbp_core_nan_cells: usize,
    rloop_resolve_core_nan_cells: usize,
    conflict_risk_core_nan_cells: usize,
    nmd_core_nan_cells: usize,
    sos_nan_cells: usize,
    rlr_nan_cells: usize,
    sii_nan_cells: usize,
    panel_coverage: Vec<SplicingInstabilityPanelCoverageJson>,
}

#[derive(Serialize)]
struct SplicingInstabilitySummaryJson {
    panel_version: &'static str,
    min_genes_per_panel_cell: usize,
    conflict_panel_enabled: bool,
    nmd_panel_enabled: bool,
    thresholds: SplicingInstabilityThresholdsJson,
    global_stats: SplicingInstabilityGlobalStatsJson,
    z_reference: SplicingInstabilityZRefJson,
    cluster_stats: Vec<SplicingInstabilityClusterStatsJson>,
    top_clusters: SplicingInstabilityTopClusterJson,
    missingness: SplicingInstabilityMissingnessJson,
}

#[derive(Serialize)]
struct PipelineStepJson {
    tool: PipelineToolJson,
    artifacts: PipelineArtifactsJson,
    cell_metrics: PipelineCellMetricsJson,
    regimes: [&'static str; 6],
}

#[derive(Serialize)]
struct PipelineToolJson {
    name: &'static str,
    stage: &'static str,
    version: &'static str,
    contract_version: &'static str,
    /// "experimental": regime/fidelity/confidence are unvalidated composites.
    signature_status: &'static str,
}

#[derive(Serialize)]
struct PipelineArtifactsJson {
    summary: &'static str,
    primary_metrics: &'static str,
    panels: &'static str,
}

#[derive(Serialize)]
struct PipelineCellMetricsJson {
    file: &'static str,
    id_column: &'static str,
    regime_column: &'static str,
    confidence_column: &'static str,
    flag_column: &'static str,
}

pub fn pipeline_out_dir(out_root: &Path) -> std::path::PathBuf {
    out_root.join(PIPELINE_DIR)
}

#[allow(clippy::too_many_arguments)]
pub fn write_pipeline_contract(
    out_dir: &Path,
    input: &InputDescriptor,
    matrix: &dyn ExpressionMatrix,
    missplicing: &MissplicingMetrics,
    imbalance: &SpliceosomeImbalanceMetrics,
    sis: &SpliceIntegrityMetrics,
    coupling: Option<&CouplingStressMetrics>,
    splicing_instability: &SplicingInstabilityMetrics,
    unspliced: Option<&UnsplicedMetrics>,
    intron_retention: Option<&IntronRetentionMetrics>,
    strata: &Strata,
    catalog: &GenesetCatalog,
) -> Result<(), InputError> {
    info!("pipeline contract: building rows");
    let rows = build_rows(matrix, missplicing, imbalance, sis, coupling)?;
    info!("pipeline contract: writing spliceqc.tsv");
    write_spliceqc_tsv(&out_dir.join("spliceqc.tsv"), &rows)?;
    info!("pipeline contract: writing panels_report.tsv");
    write_panels_report_tsv(&out_dir.join("panels_report.tsv"), matrix, catalog)?;
    info!("pipeline contract: writing summary.json");
    write_summary_json(
        &out_dir.join("summary.json"),
        input,
        &rows,
        splicing_instability,
        unspliced,
        intron_retention,
        strata,
    )?;
    info!("pipeline contract: writing pipeline_step.json");
    write_pipeline_step_json(&out_dir.join("pipeline_step.json"))?;
    Ok(())
}

fn build_rows(
    matrix: &dyn ExpressionMatrix,
    missplicing: &MissplicingMetrics,
    imbalance: &SpliceosomeImbalanceMetrics,
    sis: &SpliceIntegrityMetrics,
    coupling: Option<&CouplingStressMetrics>,
) -> Result<Vec<PipelineCellRow>, InputError> {
    let n_cells = matrix.n_cells();
    if missplicing.b_u12.len() != n_cells
        || missplicing.burden.len() != n_cells
        || missplicing.burden_star.len() != n_cells
        || imbalance.imbalance.len() != n_cells
        || imbalance.axis_u2_u1.len() != n_cells
        || sis.sis.len() != n_cells
    {
        return Err(InputError::LengthMismatch(
            "pipeline contract metric length mismatch".to_string(),
        ));
    }
    if let Some(c) = coupling
        && c.coupling_stress.len() != n_cells
    {
        return Err(InputError::LengthMismatch(
            "pipeline contract coupling length mismatch".to_string(),
        ));
    }

    let rows = (0..n_cells)
        .into_par_iter()
        .map(|cell| {
            let barcode = matrix.cell_name(cell).to_string();
            let libsize = matrix.libsize(cell);
            let nnz = matrix.nnz_cell(cell);

            // Every contract metric is a single monotone map of its source
            // metric onto [0, 1]. Non-finite sources stay non-finite: missing
            // data is written as an empty field and flagged, never coerced to 0.
            let splice_fidelity_index = unit_clamp(sis.sis[cell]);
            let intron_retention_rate = saturate(missplicing.b_u12[cell]);
            let exon_skipping_rate = saturate(imbalance.axis_u2_u1[cell].abs());
            let alt_splice_burden = rational_saturate(missplicing.burden[cell]);
            let splice_junction_noise = unit_clamp(missplicing.burden_star[cell]);
            let stress_splicing_index = coupling
                .map(|c| sigmoid(c.coupling_stress[cell]))
                .unwrap_or_else(|| saturate(imbalance.imbalance[cell]));

            let confidence = confidence_score(
                splice_fidelity_index,
                intron_retention_rate,
                alt_splice_burden,
                splice_junction_noise,
            );

            let regime = classify_regime(
                splice_fidelity_index,
                intron_retention_rate,
                exon_skipping_rate,
                alt_splice_burden,
                splice_junction_noise,
                stress_splicing_index,
            );

            let missing_metrics = [
                splice_fidelity_index,
                intron_retention_rate,
                exon_skipping_rate,
                alt_splice_burden,
                splice_junction_noise,
                stress_splicing_index,
            ]
            .iter()
            .any(|v| !v.is_finite());

            let flags = build_flags(confidence, nnz, missing_metrics);

            PipelineCellRow {
                barcode,
                sample: "unknown".to_string(),
                condition: "unknown".to_string(),
                species: "unknown".to_string(),
                libsize,
                nnz,
                expressed_genes: nnz,
                splice_fidelity_index,
                intron_retention_rate,
                exon_skipping_rate,
                alt_splice_burden,
                splice_junction_noise,
                stress_splicing_index,
                regime,
                flags,
                confidence,
            }
        })
        .collect::<Vec<_>>();
    Ok(rows)
}

fn write_spliceqc_tsv(path: &Path, rows: &[PipelineCellRow]) -> Result<(), InputError> {
    // Sort indices instead of cloning the whole row vec (was n_cells × PipelineCellRow copy).
    let mut order: Vec<usize> = (0..rows.len()).collect();
    order.sort_by(|&a, &b| rows[a].barcode.cmp(&rows[b].barcode));

    let file = File::create(path).map_err(|e| InputError::io(path, e))?;
    let mut w = BufWriter::with_capacity(1 << 16, file);
    writeln!(w, "{}", SPLICEQC_HEADER).map_err(|e| InputError::io(path, e))?;

    for &i in &order {
        let row = &rows[i];
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.barcode,
            row.sample,
            row.condition,
            row.species,
            row.libsize,
            row.nnz,
            row.expressed_genes,
            fmt_metric(row.splice_fidelity_index),
            fmt_metric(row.intron_retention_rate),
            fmt_metric(row.exon_skipping_rate),
            fmt_metric(row.alt_splice_burden),
            fmt_metric(row.splice_junction_noise),
            fmt_metric(row.stress_splicing_index),
            row.regime,
            row.flags,
            fmt_metric(row.confidence)
        )
        .map_err(|e| InputError::io(path, e))?;
    }
    w.flush().map_err(|e| InputError::io(path, e))?;
    Ok(())
}

fn write_panels_report_tsv(
    path: &Path,
    matrix: &dyn ExpressionMatrix,
    catalog: &GenesetCatalog,
) -> Result<(), InputError> {
    let rows = catalog
        .genesets
        .par_iter()
        .map(|geneset| {
            let (coverage_median, coverage_p10, sum_median, sum_p90, sum_p99) =
                panel_quantiles(matrix, geneset);
            let missing = if geneset.missing.is_empty() {
                String::new()
            } else {
                let mut missing = geneset.missing.clone();
                missing.sort();
                missing.join(",")
            };
            format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\n",
                geneset.id,
                geneset.id,
                geneset.axis,
                geneset.gene_ids.len() + geneset.missing.len(),
                geneset.gene_ids.len(),
                missing,
                coverage_median,
                coverage_p10,
                sum_median,
                sum_p90,
                sum_p99,
            )
        })
        .collect::<Vec<_>>();

    let file = File::create(path).map_err(|e| InputError::io(path, e))?;
    let mut w = BufWriter::with_capacity(1 << 16, file);
    w.write_all(b"panel_id\tpanel_name\tpanel_group\tpanel_size_defined\tpanel_size_mappable\tmissing_genes\tcoverage_median\tcoverage_p10\tsum_median\tsum_p90\tsum_p99\n")
        .map_err(|e| InputError::io(path, e))?;
    for row in rows {
        w.write_all(row.as_bytes()).map_err(|e| InputError::io(path, e))?;
    }
    w.flush().map_err(|e| InputError::io(path, e))?;
    Ok(())
}

fn panel_quantiles(matrix: &dyn ExpressionMatrix, geneset: &Geneset) -> (f64, f64, f64, f64, f64) {
    let panel_len = geneset.gene_ids.len();
    if panel_len == 0 {
        return (0.0, 0.0, 0.0, 0.0, 0.0);
    }

    let n_cells = matrix.n_cells();
    let mut coverage = Vec::with_capacity(n_cells);
    let mut sums = Vec::with_capacity(n_cells);
    let panel_len_f = panel_len as f64;
    for cell in 0..n_cells {
        // Sparse merge: O(panel + nnz_cell) instead of panel × log(nnz) bin searches.
        let (sum, detected) = matrix.panel_count_sum_and_detected(&geneset.gene_ids, cell);
        coverage.push(detected as f64 / panel_len_f);
        sums.push(sum as f64);
    }

    (
        quantile_f64(&mut coverage, 0.5),
        quantile_f64(&mut coverage, 0.1),
        quantile_f64(&mut sums, 0.5),
        quantile_f64(&mut sums, 0.9),
        quantile_f64(&mut sums, 0.99),
    )
}

fn write_summary_json(
    path: &Path,
    input: &InputDescriptor,
    rows: &[PipelineCellRow],
    splicing_instability: &SplicingInstabilityMetrics,
    unspliced: Option<&UnsplicedMetrics>,
    intron_retention: Option<&IntronRetentionMetrics>,
    strata: &Strata,
) -> Result<(), InputError> {
    let mut fidelity = rows
        .iter()
        .map(|r| r.splice_fidelity_index)
        .filter(|v| v.is_finite())
        .collect::<Vec<_>>();
    let mut stress = rows
        .iter()
        .map(|r| r.stress_splicing_index)
        .filter(|v| v.is_finite())
        .collect::<Vec<_>>();

    let mut counts: BTreeMap<&'static str, u64> = BTreeMap::new();
    for regime in REGIMES {
        counts.insert(regime, 0);
    }
    for row in rows {
        if let Some(entry) = counts.get_mut(row.regime) {
            *entry += 1;
        }
    }

    let n = rows.len().max(1) as f64;
    let mut fractions = BTreeMap::new();
    for regime in REGIMES {
        let value = *counts.get(regime).unwrap_or(&0) as f64 / n;
        fractions.insert(regime, value);
    }

    let low_conf = rows.iter().filter(|r| r.confidence < 0.5).count() as f64 / n;
    let high_noise = rows
        .iter()
        .filter(|r| r.splice_junction_noise > 0.7)
        .count() as f64
        / n;

    let summary = SummaryJson {
        tool: ToolJson {
            name: "kira-spliceqc",
            version: env!("CARGO_PKG_VERSION"),
            simd: crate::simd::backend(),
            experimental_signatures: true,
        },
        input: InputJson {
            n_cells: input.n_cells,
            species: "unknown",
            levels: if unspliced.is_some() { vec!["L0", "L1"] } else { vec!["L0"] },
        },
        distributions: DistributionsJson {
            splice_fidelity_index: DistributionJson {
                median: quantile_f64(&mut fidelity, 0.5),
                p90: quantile_f64(&mut fidelity, 0.9),
                p99: quantile_f64(&mut fidelity, 0.99),
            },
            stress_splicing_index: DistributionJson {
                median: quantile_f64(&mut stress, 0.5),
                p90: quantile_f64(&mut stress, 0.9),
                p99: quantile_f64(&mut stress, 0.99),
            },
        },
        regimes: RegimesJson { counts, fractions },
        qc: QcJson {
            low_confidence_fraction: low_conf,
            high_splice_noise_fraction: high_noise,
        },
        splicing_instability: build_splicing_instability_summary(splicing_instability),
        unspliced: unspliced.map(|u| {
            let mut defined: Vec<f64> = u
                .unspliced_fraction
                .iter()
                .filter(|v| v.is_finite())
                .map(|v| *v as f64)
                .collect();
            let opt = |v: f64| if v.is_finite() { Some(v) } else { None };
            UnsplicedSummaryJson {
                source: u.source.clone(),
                min_layer_umis: u.min_layer_umis,
                n_defined_cells: defined.len(),
                cells_without_layers: u.cells_without_layers,
                median: opt(quantile_f64(&mut defined, 0.5)),
                p10: opt(quantile_f64(&mut defined, 0.1)),
                p90: opt(quantile_f64(&mut defined, 0.9)),
                nuclear_fraction_flag_fraction: u.nuclear_fraction_flag.iter().filter(|f| **f).count()
                    as f64
                    / u.nuclear_fraction_flag.len().max(1) as f64,
                strata: u
                    .reference
                    .iter()
                    .map(|s| StratumStatJson {
                        name: s.name.clone(),
                        n_cells: s.n_cells,
                        n_defined: s.n_defined,
                        median: opt_f32(s.median),
                        mad: opt_f32(s.mad),
                    })
                    .collect(),
            }
        }),
        intron_retention: intron_retention.map(|m| {
            let mut defined: Vec<f64> = m
                .intron_retention_index
                .iter()
                .filter(|v| v.is_finite())
                .map(|v| *v as f64)
                .collect();
            let opt = |v: f64| if v.is_finite() { Some(v) } else { None };
            IntronRetentionSummaryJson {
                min_gene_umis: m.min_gene_umis,
                min_genes: m.min_genes,
                genes_with_reference: m.genes_with_reference,
                n_defined_cells: defined.len(),
                median: opt(quantile_f64(&mut defined, 0.5)),
                p10: opt(quantile_f64(&mut defined, 0.1)),
                p90: opt(quantile_f64(&mut defined, 0.9)),
                high_flag_fraction: m.intron_retention_high.iter().filter(|f| **f).count() as f64
                    / m.intron_retention_high.len().max(1) as f64,
                strata: m
                    .reference
                    .iter()
                    .map(|s| StratumStatJson {
                        name: s.name.clone(),
                        n_cells: s.n_cells,
                        n_defined: s.n_defined,
                        median: opt_f32(s.median),
                        mad: opt_f32(s.mad),
                    })
                    .collect(),
            }
        }),
        reference: ReferenceJson {
            mode: strata.mode.as_str(),
            column: strata.column.clone(),
            n_strata: strata.n_strata(),
            min_stratum_cells: MIN_STRATUM_CELLS,
            folded_cells: strata.folded_cells,
            strata: strata
                .names
                .iter()
                .zip(strata.sizes())
                .map(|(name, n_cells)| StratumSizeJson {
                    name: name.clone(),
                    n_cells,
                })
                .collect(),
        },
    };

    let file = File::create(path).map_err(|e| InputError::io(path, e))?;
    let mut w = BufWriter::with_capacity(1 << 16, file);
    serde_json::to_writer(&mut w, &summary)
        .map_err(|e| InputError::OutputSerialization(e.to_string()))?;
    w.flush().map_err(|e| InputError::io(path, e))?;
    Ok(())
}

fn build_splicing_instability_summary(
    metrics: &SplicingInstabilityMetrics,
) -> SplicingInstabilitySummaryJson {
    let mut by_sii = metrics
        .cluster_stats
        .iter()
        .filter(|c| c.sii_median.is_finite())
        .collect::<Vec<_>>();
    by_sii.sort_by(|a, b| {
        b.sii_median
            .partial_cmp(&a.sii_median)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.cluster_id.cmp(&b.cluster_id))
    });

    let mut by_flag = metrics.cluster_stats.iter().collect::<Vec<_>>();
    by_flag.sort_by(|a, b| {
        b.genome_instability_splicing_flag_fraction
            .partial_cmp(&a.genome_instability_splicing_flag_fraction)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.cluster_id.cmp(&b.cluster_id))
    });

    SplicingInstabilitySummaryJson {
        panel_version: metrics.panel_version,
        min_genes_per_panel_cell: metrics.min_genes,
        conflict_panel_enabled: metrics.conflict_panel_enabled,
        nmd_panel_enabled: metrics.nmd_panel_enabled,
        thresholds: SplicingInstabilityThresholdsJson {
            splice_overload_high: SPLICE_OVERLOAD_HIGH_THRESHOLD,
            rloop_risk_high: RLOOP_RISK_HIGH_THRESHOLD,
            splicing_instability_high: SPLICING_INSTABILITY_HIGH_THRESHOLD,
        },
        global_stats: SplicingInstabilityGlobalStatsJson {
            sos_p50: opt_f32(metrics.global_stats.sos_p50),
            sos_p90: opt_f32(metrics.global_stats.sos_p90),
            rlr_p50: opt_f32(metrics.global_stats.rlr_p50),
            rlr_p90: opt_f32(metrics.global_stats.rlr_p90),
            sii_p50: opt_f32(metrics.global_stats.sii_p50),
            sii_p90: opt_f32(metrics.global_stats.sii_p90),
        },
        z_reference: SplicingInstabilityZRefJson {
            splice_core: SplicingInstabilityRobustRefJson {
                median: opt_f32(metrics.z_reference.splice_core.median),
                mad: opt_f32(metrics.z_reference.splice_core.mad),
            },
            rbp_core: SplicingInstabilityRobustRefJson {
                median: opt_f32(metrics.z_reference.rbp_core.median),
                mad: opt_f32(metrics.z_reference.rbp_core.mad),
            },
            rloop_resolve_core: SplicingInstabilityRobustRefJson {
                median: opt_f32(metrics.z_reference.rloop_resolve_core.median),
                mad: opt_f32(metrics.z_reference.rloop_resolve_core.mad),
            },
            conflict_risk_core: metrics.z_reference.conflict_risk_core.as_ref().map(|r| {
                SplicingInstabilityRobustRefJson {
                    median: opt_f32(r.median),
                    mad: opt_f32(r.mad),
                }
            }),
            nmd_core: metrics.z_reference.nmd_core.as_ref().map(|r| {
                SplicingInstabilityRobustRefJson {
                    median: opt_f32(r.median),
                    mad: opt_f32(r.mad),
                }
            }),
        },
        cluster_stats: metrics
            .cluster_stats
            .iter()
            .map(|cluster| SplicingInstabilityClusterStatsJson {
                cluster_id: cluster.cluster_id.clone(),
                n_cells: cluster.n_cells,
                sos_median: opt_f32(cluster.sos_median),
                sos_p10: opt_f32(cluster.sos_p10),
                sos_p90: opt_f32(cluster.sos_p90),
                rlr_median: opt_f32(cluster.rlr_median),
                rlr_p10: opt_f32(cluster.rlr_p10),
                rlr_p90: opt_f32(cluster.rlr_p90),
                sii_median: opt_f32(cluster.sii_median),
                sii_p10: opt_f32(cluster.sii_p10),
                sii_p90: opt_f32(cluster.sii_p90),
                splice_overload_high_fraction: cluster.splice_overload_high_fraction,
                rloop_risk_high_fraction: cluster.rloop_risk_high_fraction,
                splicing_instability_high_fraction: cluster.splicing_instability_high_fraction,
                genome_instability_splicing_flag_fraction: cluster
                    .genome_instability_splicing_flag_fraction,
            })
            .collect(),
        top_clusters: SplicingInstabilityTopClusterJson {
            by_sii_median: by_sii
                .into_iter()
                .take(5)
                .map(|cluster| cluster.cluster_id.clone())
                .collect(),
            by_genome_instability_splicing_flag_fraction: by_flag
                .into_iter()
                .take(5)
                .map(|cluster| cluster.cluster_id.clone())
                .collect(),
        },
        missingness: SplicingInstabilityMissingnessJson {
            splice_core_nan_cells: metrics.missingness.splice_core_nan_cells,
            rbp_core_nan_cells: metrics.missingness.rbp_core_nan_cells,
            rloop_resolve_core_nan_cells: metrics.missingness.rloop_resolve_core_nan_cells,
            conflict_risk_core_nan_cells: metrics.missingness.conflict_risk_core_nan_cells,
            nmd_core_nan_cells: metrics.missingness.nmd_core_nan_cells,
            sos_nan_cells: metrics.missingness.sos_nan_cells,
            rlr_nan_cells: metrics.missingness.rlr_nan_cells,
            sii_nan_cells: metrics.missingness.sii_nan_cells,
            panel_coverage: metrics
                .missingness
                .panel_coverage
                .iter()
                .map(|panel| SplicingInstabilityPanelCoverageJson {
                    panel_name: panel.panel_name,
                    genes_defined: panel.genes_defined,
                    genes_mapped: panel.genes_mapped,
                })
                .collect(),
        },
    }
}

fn write_pipeline_step_json(path: &Path) -> Result<(), InputError> {
    let step = PipelineStepJson {
        tool: PipelineToolJson {
            name: "kira-spliceqc",
            stage: "splicing",
            version: env!("CARGO_PKG_VERSION"),
            contract_version: PIPELINE_CONTRACT_VERSION,
            signature_status: "experimental",
        },
        artifacts: PipelineArtifactsJson {
            summary: "summary.json",
            primary_metrics: "spliceqc.tsv",
            panels: "panels_report.tsv",
        },
        cell_metrics: PipelineCellMetricsJson {
            file: "spliceqc.tsv",
            id_column: "barcode",
            regime_column: "regime",
            confidence_column: "confidence",
            flag_column: "flags",
        },
        regimes: REGIMES,
    };
    let file = File::create(path).map_err(|e| InputError::io(path, e))?;
    let mut w = BufWriter::with_capacity(1 << 16, file);
    serde_json::to_writer(&mut w, &step)
        .map_err(|e| InputError::OutputSerialization(e.to_string()))?;
    w.flush().map_err(|e| InputError::io(path, e))?;
    Ok(())
}

pub fn spliceqc_header() -> &'static str {
    SPLICEQC_HEADER
}

/// Identity on `[0, 1]` for metrics that are unit-bounded by construction
/// (SIS, burden*). Values outside the range are clamped; NaN is preserved.
fn unit_clamp(value: f32) -> f64 {
    if !value.is_finite() {
        return f64::NAN;
    }
    (value as f64).clamp(0.0, 1.0)
}

/// `1 - exp(-x)` for non-negative, unbounded burden-like inputs. Strictly
/// monotone on `[0, inf)`, maps 0 -> 0. Negative inputs are treated as 0.
fn saturate(value: f32) -> f64 {
    if !value.is_finite() {
        return f64::NAN;
    }
    let x = (value as f64).max(0.0);
    1.0 - (-x).exp()
}

/// `x / (1 + x)` for non-negative inputs. Strictly monotone on `[0, inf)`.
/// Used where a second, distinct saturation of the same burden is part of
/// the deprecated contract (see METRICS.md).
fn rational_saturate(value: f32) -> f64 {
    if !value.is_finite() {
        return f64::NAN;
    }
    let x = (value as f64).max(0.0);
    x / (1.0 + x)
}

/// Logistic map for signed inputs such as `coupling_stress`.
fn sigmoid(value: f32) -> f64 {
    if !value.is_finite() {
        return f64::NAN;
    }
    let x = value as f64;
    1.0 / (1.0 + (-x).exp())
}

/// Six-decimal fixed formatting; non-finite values become an empty field.
fn fmt_metric(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.6}")
    } else {
        String::new()
    }
}

fn confidence_score(fidelity: f64, intron_retention: f64, alt_splice: f64, noise: f64) -> f64 {
    if !(fidelity.is_finite()
        && intron_retention.is_finite()
        && alt_splice.is_finite()
        && noise.is_finite())
    {
        return f64::NAN;
    }
    let penalty = 0.4 * intron_retention + 0.3 * alt_splice + 0.3 * noise;
    (0.65 * fidelity + 0.35 * (1.0 - penalty)).clamp(0.0, 1.0)
}

fn classify_regime(
    fidelity: f64,
    intron_retention: f64,
    exon_skipping: f64,
    alt_splice: f64,
    noise: f64,
    stress: f64,
) -> &'static str {
    if !(fidelity.is_finite()
        && intron_retention.is_finite()
        && exon_skipping.is_finite()
        && alt_splice.is_finite()
        && noise.is_finite()
        && stress.is_finite())
    {
        return "Unclassified";
    }
    if fidelity < 0.2 && stress > 0.8 {
        "SplicingCollapse"
    } else if noise > 0.75 {
        "SpliceNoiseDominant"
    } else if stress > 0.65 {
        "StressInducedSplicing"
    } else if alt_splice > 0.45 || exon_skipping > 0.45 || intron_retention > 0.45 {
        "RegulatedAlternativeSplicing"
    } else if fidelity > 0.75 && noise < 0.35 && stress < 0.4 {
        "HighFidelitySplicing"
    } else {
        "Unclassified"
    }
}

fn build_flags(confidence: f64, nnz: u64, missing_metrics: bool) -> String {
    let mut flags = Vec::new();
    if confidence.is_finite() && confidence < 0.5 {
        flags.push("LOW_CONFIDENCE");
    }
    if nnz < 50 {
        flags.push("LOW_SPLICE_SIGNAL");
    }
    if missing_metrics {
        flags.push("MISSING_METRICS");
    }
    flags.join(",")
}

fn opt_f32(value: f32) -> Option<f32> {
    if value.is_finite() { Some(value) } else { None }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_monotone(name: &str, f: fn(f32) -> f64, lo: f32, hi: f32) {
        let steps = 4000;
        let mut prev = f(lo);
        assert!(prev.is_finite(), "{name}: non-finite at lower bound");
        for i in 1..=steps {
            let x = lo + (hi - lo) * (i as f32 / steps as f32);
            let y = f(x);
            assert!(y.is_finite(), "{name}: non-finite at {x}");
            assert!(
                y >= prev,
                "{name}: not monotone between {} and {x} ({prev} > {y})",
                lo + (hi - lo) * ((i - 1) as f32 / steps as f32)
            );
            assert!((0.0..=1.0).contains(&y), "{name}: {y} outside [0, 1]");
            prev = y;
        }
    }

    #[test]
    fn contract_transforms_are_monotone_and_unit_bounded() {
        assert_monotone("unit_clamp", unit_clamp, -2.0, 3.0);
        assert_monotone("saturate", saturate, -2.0, 20.0);
        assert_monotone("rational_saturate", rational_saturate, -2.0, 20.0);
        assert_monotone("sigmoid", sigmoid, -12.0, 12.0);
    }

    #[test]
    fn contract_transforms_preserve_nan() {
        for f in [unit_clamp, saturate, rational_saturate, sigmoid] {
            assert!(f(f32::NAN).is_nan());
            assert!(f(f32::INFINITY).is_nan());
        }
        assert!(confidence_score(f64::NAN, 0.1, 0.1, 0.1).is_nan());
        assert_eq!(
            classify_regime(f64::NAN, 0.1, 0.1, 0.1, 0.1, 0.1),
            "Unclassified"
        );
        assert_eq!(fmt_metric(f64::NAN), "");
        assert_eq!(fmt_metric(0.5), "0.500000");
    }

    #[test]
    fn missing_metrics_flag_is_set_without_low_confidence() {
        let flags = build_flags(f64::NAN, 1000, true);
        assert_eq!(flags, "MISSING_METRICS");
    }
}
