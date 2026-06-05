use std::collections::HashSet;
use std::time::Instant;

use rayon::prelude::*;
use tracing::{debug, info, warn};

use crate::expression::ExpressionMatrix;
use crate::genesets::catalog::default_catalog_path;
use crate::genesets::{GenesetCatalog, load_catalog};
use crate::input::error::InputError;
use crate::model::isoform_dispersion::IsoformDispersionMetrics;
use crate::stats::robust::{mad, median};

const EPS: f32 = 1e-12;

const REGULATOR_SET_IDS: &[&str] = &[
    "SRSF_SR", "HNRNP", "U1_CORE", "U2_CORE", "SF3B_AXIS", "MINOR_U12",
];

pub fn run_stage3(matrix: &dyn ExpressionMatrix) -> Result<IsoformDispersionMetrics, InputError> {
    let catalog_path = default_catalog_path();
    let catalog = load_catalog(&catalog_path, matrix)?;
    compute(matrix, &catalog)
}

pub fn compute(
    matrix: &dyn ExpressionMatrix,
    catalog: &GenesetCatalog,
) -> Result<IsoformDispersionMetrics, InputError> {
    let regulator_genes = build_regulator_union(catalog)?;
    let n_reg = regulator_genes.len();
    info!(regulators = n_reg, "regulator union size");

    let n_cells = matrix.n_cells();
    let ln_r = (n_reg as f32).ln();

    let start = Instant::now();

    // One sparse pass per cell collects cp10k(reg, cell) values, then computes
    // sum / entropy / dispersion in the same pass. Was: two passes × bin-search.
    let mut entropy = vec![f32::NAN; n_cells];
    let mut dispersion = vec![f32::NAN; n_cells];

    entropy
        .par_iter_mut()
        .zip(dispersion.par_iter_mut())
        .enumerate()
        .for_each_init(
            || Vec::<f32>::with_capacity(n_reg),
            |scratch, (cell, (ent_out, disp_out))| {
                scratch.clear();
                // cp10k = log_cp10k.exp() - 1 — but we want cp10k itself. Use the
                // matrix's gather (log_cp10k) and back-transform with expm1 to
                // recover cp10k without a second binary search.
                matrix.gather_panel_log_cp10k(&regulator_genes, cell, scratch);

                let mut sum = 0.0_f32;
                for v in scratch.iter() {
                    sum += v.exp() - 1.0;
                }
                if sum <= 0.0 {
                    *ent_out = f32::NAN;
                    *disp_out = f32::NAN;
                    return;
                }

                let mut h = 0.0_f32;
                let mut q = 0.0_f32;
                let inv_sum_eps = 1.0 / (sum + EPS);
                for v in scratch.iter() {
                    let y = v.exp() - 1.0;
                    let p = y * inv_sum_eps;
                    h -= p * (p + EPS).ln();
                    q += p * p;
                }
                *ent_out = if ln_r > 0.0 { h / ln_r } else { f32::NAN };
                *disp_out = (1.0 / (q + EPS)) / n_reg as f32;
            },
        );

    let zero_count = entropy.iter().filter(|v| v.is_nan()).count();
    if zero_count > 0 {
        warn!(
            zero_cells = zero_count,
            "cells with zero regulator expression or undefined entropy"
        );
    }

    let med = median(&entropy);
    let mad_val = mad(&entropy, med) * 1.4826 + EPS;
    let z_entropy: Vec<f32> = entropy
        .iter()
        .map(|&value| {
            if value.is_finite() && med.is_finite() && mad_val.is_finite() {
                (value - med) / mad_val
            } else {
                f32::NAN
            }
        })
        .collect();

    debug!(
        elapsed_ms = start.elapsed().as_millis(),
        "isoform dispersion computed"
    );

    Ok(IsoformDispersionMetrics {
        entropy,
        dispersion,
        z_entropy,
    })
}

fn build_regulator_union(catalog: &GenesetCatalog) -> Result<Vec<u32>, InputError> {
    let mut ids: HashSet<u32> = HashSet::new();
    for geneset in &catalog.genesets {
        if REGULATOR_SET_IDS.contains(&geneset.id.as_str()) {
            for &gene_id in &geneset.gene_ids {
                ids.insert(gene_id);
            }
        }
    }
    if ids.is_empty() {
        return Err(InputError::EmptyRegulatorSet);
    }
    let mut ids_vec: Vec<u32> = ids.into_iter().collect();
    ids_vec.sort_unstable();
    Ok(ids_vec)
}
