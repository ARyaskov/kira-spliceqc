//! Provenance block written to `summary.json` and `cells.json`: everything
//! needed to reproduce or audit a run without the log.

use std::path::Path;

use serde::Serialize;

use crate::cli::config::{RunConfig, RunMode};
use crate::genesets::controls::CONTROLS_PER_GENE;
use crate::input::shared_cache::crc64_ecma;
use crate::metrics::intron_retention::{
    MIN_CELLS_PER_GENE, MIN_GENE_UMIS, MIN_GENES as IRI_MIN_GENES, PRIOR_STRENGTH, WEIGHT_CAP_UMIS,
};
use crate::metrics::junctions::{CRYPTIC_MAX, CRYPTIC_MIN, MIN_JUNCTION_UMIS, MIN_RATIO_UMIS};
use crate::metrics::unspliced::MIN_LAYER_UMIS;
use crate::model::cell_cycle::CellCycleMetrics;
use crate::model::cell_qc::CellQc;
use crate::model::intron_retention::IntronRetentionMetrics;
use crate::model::junctions::JunctionMetrics;
use crate::model::sis::SpliceIntegrityMetrics;
use crate::model::splicing_instability::SplicingInstabilityMetrics;
use crate::model::unspliced::UnsplicedMetrics;
use crate::reference::{DEVIATION_THRESHOLD, FLAG_FDR, MAX_DEPTH_BINS, MIN_STRATUM_CELLS, Strata};

#[derive(Debug, Clone, Serialize)]
pub struct Provenance {
    pub tool: ToolInfo,
    pub command: CommandInfo,
    pub geneset_catalog: FileInfo,
    pub instability_panel_version: &'static str,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reference_file: Option<FileInfo>,
    pub input_levels: Vec<&'static str>,
    pub reference: ReferenceInfo,
    pub parameters: Parameters,
    /// Cells with an undefined value per metric.
    pub undefined_cells: UndefinedCells,
}

#[derive(Debug, Clone, Serialize)]
pub struct ToolInfo {
    pub name: &'static str,
    pub version: &'static str,
    pub simd: &'static str,
}

#[derive(Debug, Clone, Serialize)]
pub struct CommandInfo {
    pub input: String,
    pub run_mode: &'static str,
    pub extended: bool,
    pub experimental_signatures: bool,
    pub layers: Option<String>,
    pub metadata: Option<String>,
    pub stratify_by: Option<String>,
    pub reference: Option<String>,
    pub threads: Option<usize>,
}

#[derive(Debug, Clone, Serialize)]
pub struct FileInfo {
    pub source: String,
    /// CRC-64/ECMA of the file bytes (hex).
    pub crc64: String,
    pub bytes: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct ReferenceInfo {
    pub mode: &'static str,
    pub column: Option<String>,
    pub n_strata: usize,
    pub folded_cells: usize,
    /// Cells excluded from norm computation (low depth or doublet).
    pub excluded_cells: usize,
    pub doublet_column: Option<String>,
    /// The external reference also provided expression-signature norms.
    pub external_expression_norms: bool,
}

#[derive(Debug, Clone, Serialize)]
pub struct Parameters {
    pub min_junction_umis: u64,
    pub min_ratio_umis: u64,
    pub cryptic_window_nt: [u64; 2],
    pub min_counts: u64,
    pub min_genes: u64,
    pub controls_per_gene: usize,
    pub min_stratum_cells: usize,
    pub max_depth_bins: usize,
    pub deviation_threshold: f32,
    pub flag_fdr: f64,
    pub min_layer_umis: u64,
    pub iri_min_gene_umis: u32,
    pub iri_min_genes: usize,
    pub iri_min_cells_per_gene: usize,
    pub iri_prior_strength: f64,
    pub iri_weight_cap_umis: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct UndefinedCells {
    pub sis: usize,
    pub sos: usize,
    pub unspliced_fraction: Option<usize>,
    pub intron_retention_index: Option<usize>,
    pub junction_metrics: Option<usize>,
    pub cell_cycle_phase: usize,
}

/// `["L0"]`, `["L0", "L1"]`, `["L0", "L2"]` or `["L0", "L1", "L2"]`.
pub fn input_levels(has_layers: bool, has_junctions: bool) -> Vec<&'static str> {
    let mut v = vec!["L0"];
    if has_layers {
        v.push("L1");
    }
    if has_junctions {
        v.push("L2");
    }
    v
}

impl FileInfo {
    pub fn of_bytes(source: impl Into<String>, bytes: &[u8]) -> Self {
        Self {
            source: source.into(),
            crc64: format!("{:016x}", crc64_ecma(bytes)),
            bytes: bytes.len(),
        }
    }

    pub fn of_path(path: &Path) -> Option<Self> {
        std::fs::read(path)
            .ok()
            .map(|bytes| Self::of_bytes(path.display().to_string(), &bytes))
    }
}

#[allow(clippy::too_many_arguments)]
pub fn build(
    config: &RunConfig,
    geneset_catalog: FileInfo,
    reference_file: Option<FileInfo>,
    external_expression_norms: bool,
    has_layers: bool,
    junctions: Option<&JunctionMetrics>,
    strata: &Strata,
    cell_qc: &CellQc,
    sis: &SpliceIntegrityMetrics,
    instability: &SplicingInstabilityMetrics,
    unspliced: Option<&UnsplicedMetrics>,
    intron_retention: Option<&IntronRetentionMetrics>,
    cell_cycle: &CellCycleMetrics,
) -> Provenance {
    let count_nan = |v: &[f32]| v.iter().filter(|x| !x.is_finite()).count();
    Provenance {
        tool: ToolInfo {
            name: "kira-spliceqc",
            version: env!("CARGO_PKG_VERSION"),
            simd: crate::simd::backend(),
        },
        command: CommandInfo {
            input: config.input.display().to_string(),
            run_mode: match config.run_mode {
                RunMode::Standalone => "standalone",
                RunMode::Pipeline => "pipeline",
            },
            extended: config.extended,
            experimental_signatures: config.experimental_signatures
                || config.run_mode == RunMode::Pipeline,
            layers: config.layers.as_ref().map(|p| p.display().to_string()),
            metadata: config.metadata.as_ref().map(|p| p.display().to_string()),
            stratify_by: config.stratify_by.clone(),
            reference: config.reference.as_ref().map(|p| p.display().to_string()),
            threads: config.threads,
        },
        geneset_catalog,
        instability_panel_version: instability.panel_version,
        reference_file,
        input_levels: input_levels(has_layers, junctions.is_some()),
        reference: ReferenceInfo {
            mode: strata.mode.as_str(),
            column: strata.column.clone(),
            n_strata: strata.n_strata(),
            folded_cells: strata.folded_cells,
            excluded_cells: strata.n_excluded(),
            doublet_column: cell_qc.doublet_column.clone(),
            external_expression_norms,
        },
        parameters: Parameters {
            min_junction_umis: MIN_JUNCTION_UMIS,
            min_ratio_umis: MIN_RATIO_UMIS,
            cryptic_window_nt: [CRYPTIC_MIN, CRYPTIC_MAX],
            min_counts: cell_qc.min_counts,
            min_genes: cell_qc.min_genes,
            controls_per_gene: CONTROLS_PER_GENE,
            min_stratum_cells: MIN_STRATUM_CELLS,
            max_depth_bins: MAX_DEPTH_BINS,
            deviation_threshold: DEVIATION_THRESHOLD,
            flag_fdr: FLAG_FDR,
            min_layer_umis: MIN_LAYER_UMIS,
            iri_min_gene_umis: MIN_GENE_UMIS,
            iri_min_genes: IRI_MIN_GENES,
            iri_min_cells_per_gene: MIN_CELLS_PER_GENE,
            iri_prior_strength: PRIOR_STRENGTH,
            iri_weight_cap_umis: WEIGHT_CAP_UMIS,
        },
        undefined_cells: UndefinedCells {
            sis: count_nan(&sis.sis),
            sos: count_nan(&instability.sos),
            unspliced_fraction: unspliced.map(|u| u.undefined_cells),
            intron_retention_index: intron_retention.map(|m| m.undefined_cells),
            junction_metrics: junctions.map(|m| m.undefined_cells),
            cell_cycle_phase: cell_cycle
                .phase
                .iter()
                .filter(|p| **p == crate::model::cell_cycle::CellCyclePhase::Unknown)
                .count(),
        },
    }
}
