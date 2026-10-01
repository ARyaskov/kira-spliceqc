//! External reference files (`ref.json`).
//!
//! A reference is built from a control dataset (`kira-spliceqc reference
//! build`) and applied to another run with `--reference`. It carries, per
//! stratum, the Tier A norms: the logit median and overdispersion of the
//! unspliced fraction, the median and overdispersion of the intron retention
//! index, and the pooled per-gene unspliced ratios. Cells of the target run
//! are assigned to reference strata by the same metadata column; cells whose
//! value has no reference stratum fall back to the reference's `global`
//! stratum. Expression signatures stay dataset-relative (their depth-binned
//! standardization has no closed-form external norm yet); `summary.json`
//! records which metrics used the external norms.

use std::collections::BTreeMap;
use std::path::Path;

use serde::{Deserialize, Serialize};
use tracing::{info, warn};

use crate::expression::ExpressionMatrix;
use crate::genesets::aliases::{resolve_symbol, symbol_index};
use crate::input::error::InputError;
use crate::input::metadata::CellMetadata;
use crate::model::intron_retention::IntronRetentionMetrics;
use crate::model::unspliced::UnsplicedMetrics;
use crate::reference::{ContinuousNorm, GLOBAL_STRATUM, ProportionNorm, ReferenceMode, Strata};

pub const REFERENCE_FORMAT: &str = "kira-spliceqc-reference";
pub const REFERENCE_VERSION: u32 = 1;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceFile {
    pub format: String,
    pub version: u32,
    pub tool_version: String,
    /// Metadata column that defines the strata (`null` = single global stratum).
    pub stratify_by: Option<String>,
    pub n_cells: usize,
    pub min_layer_umis: u64,
    /// Strata in label order; the first is always `global`.
    pub strata: Vec<ReferenceStratum>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceStratum {
    pub name: String,
    pub n_cells: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub unspliced_fraction: Option<ProportionNorm>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub intron_retention_index: Option<ContinuousNorm>,
    /// Pooled unspliced ratio per gene symbol (genes with a defined reference).
    #[serde(default)]
    pub gene_unspliced_ratio: BTreeMap<String, f32>,
}

impl ReferenceFile {
    pub fn read(path: &Path) -> Result<Self, InputError> {
        let bytes = std::fs::read(path).map_err(|e| InputError::io(path, e))?;
        let file: ReferenceFile = serde_json::from_slice(&bytes)
            .map_err(|e| InputError::InvalidReference(format!("{}: {e}", path.display())))?;
        if file.format != REFERENCE_FORMAT {
            return Err(InputError::InvalidReference(format!(
                "{}: format {:?}, expected {REFERENCE_FORMAT:?}",
                path.display(),
                file.format
            )));
        }
        if file.version != REFERENCE_VERSION {
            return Err(InputError::InvalidReference(format!(
                "{}: version {}, this build reads version {REFERENCE_VERSION}",
                path.display(),
                file.version
            )));
        }
        if file.strata.is_empty() || file.strata[0].name != GLOBAL_STRATUM {
            return Err(InputError::InvalidReference(format!(
                "{}: the first stratum must be {GLOBAL_STRATUM:?}",
                path.display()
            )));
        }
        info!(path = %path.display(), strata = file.strata.len(), "external reference loaded");
        Ok(file)
    }

    pub fn write(&self, path: &Path) -> Result<(), InputError> {
        let file = std::fs::File::create(path).map_err(|e| InputError::io(path, e))?;
        let mut w = std::io::BufWriter::new(file);
        serde_json::to_writer_pretty(&mut w, self)
            .map_err(|e| InputError::OutputSerialization(e.to_string()))?;
        use std::io::Write;
        w.flush().map_err(|e| InputError::io(path, e))
    }

    /// Assigns the cells of a target run to the reference strata by the
    /// reference's metadata column; unmatched cells go to `global`.
    pub fn assign(&self, metadata: &CellMetadata, n_cells: usize) -> Strata {
        let names: Vec<String> = self.strata.iter().map(|s| s.name.clone()).collect();
        let index: BTreeMap<&str, u32> = names
            .iter()
            .enumerate()
            .map(|(i, n)| (n.as_str(), i as u32))
            .collect();
        let (labels, folded, column) = match self.stratify_by.as_deref().and_then(|col| {
            metadata
                .resolve(&[col])
                .map(|(found, values)| (found.to_string(), values))
        }) {
            Some((found, values)) if values.len() == n_cells => {
                let mut folded = 0usize;
                let labels = values
                    .iter()
                    .map(|v| match index.get(v.as_str()) {
                        Some(&l) => l,
                        None => {
                            folded += 1;
                            0
                        }
                    })
                    .collect();
                (labels, folded, Some(found))
            }
            _ => {
                if let Some(col) = &self.stratify_by {
                    warn!(
                        column = col.as_str(),
                        "reference stratification column not in this dataset; every cell uses the global reference stratum"
                    );
                }
                (
                    vec![0u32; n_cells],
                    if self.stratify_by.is_some() {
                        n_cells
                    } else {
                        0
                    },
                    None,
                )
            }
        };
        Strata {
            mode: ReferenceMode::External,
            column,
            labels,
            names,
            folded_cells: folded,
            excluded: vec![false; n_cells],
        }
    }

    pub fn unspliced_norms(&self) -> Vec<Option<ProportionNorm>> {
        self.strata.iter().map(|s| s.unspliced_fraction).collect()
    }

    pub fn intron_retention_norms(&self) -> Vec<Option<ContinuousNorm>> {
        self.strata
            .iter()
            .map(|s| s.intron_retention_index)
            .collect()
    }

    /// Per-stratum pooled unspliced ratio per gene id of `matrix` (NaN where
    /// the reference has no ratio for the gene).
    pub fn gene_ratios_for(&self, matrix: &dyn ExpressionMatrix) -> Vec<Vec<f64>> {
        let index = symbol_index(matrix);
        self.strata
            .iter()
            .map(|s| {
                let mut p = vec![f64::NAN; matrix.n_genes()];
                for (symbol, ratio) in &s.gene_unspliced_ratio {
                    if let Some((id, _)) = resolve_symbol(&index, symbol) {
                        p[id as usize] = *ratio as f64;
                    }
                }
                p
            })
            .collect()
    }
}

/// Builds a reference from a run's strata and Tier A metrics.
pub fn build_reference(
    strata: &Strata,
    n_cells: usize,
    unspliced: Option<&UnsplicedMetrics>,
    intron_retention: Option<&IntronRetentionMetrics>,
    gene_symbols: &[String],
) -> ReferenceFile {
    let sizes = strata.sizes();
    let unspliced_norms: Vec<Option<ProportionNorm>> = match unspliced {
        Some(u) => u.norms.iter().map(|n| Some(*n)).collect(),
        None => vec![None; strata.n_strata()],
    };
    let ir_norms: Vec<Option<ContinuousNorm>> = match intron_retention {
        Some(ir) => ir.norms.iter().map(|n| Some(*n)).collect(),
        None => vec![None; strata.n_strata()],
    };
    let strata_out = strata
        .names
        .iter()
        .enumerate()
        .map(|(label, name)| {
            let gene_unspliced_ratio = intron_retention
                .map(|ir| {
                    ir.gene_reference[label]
                        .iter()
                        .enumerate()
                        .filter(|(_, p)| p.is_finite())
                        .map(|(g, p)| (gene_symbols[g].clone(), *p as f32))
                        .collect()
                })
                .unwrap_or_default();
            ReferenceStratum {
                name: name.clone(),
                n_cells: sizes[label],
                unspliced_fraction: unspliced_norms[label],
                intron_retention_index: ir_norms[label],
                gene_unspliced_ratio,
            }
        })
        .collect();
    ReferenceFile {
        format: REFERENCE_FORMAT.to_string(),
        version: REFERENCE_VERSION,
        tool_version: env!("CARGO_PKG_VERSION").to_string(),
        stratify_by: strata.column.clone(),
        n_cells,
        min_layer_umis: unspliced.map_or(0, |u| u.min_layer_umis),
        strata: strata_out,
    }
}
