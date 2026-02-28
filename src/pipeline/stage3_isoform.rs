use std::collections::HashSet;
use std::time::Instant;

use tracing::{debug, info, warn};

use crate::expression::ExpressionMatrix;
use crate::genesets::catalog::default_catalog_path;
use crate::genesets::{GenesetCatalog, load_catalog};
use crate::input::error::InputError;
use crate::model::isoform_dispersion::IsoformDispersionMetrics;
use crate::stats::robust::{mad, median};

const EPS: f32 = 1e-12;

const REGULATOR_SET_IDS: &[&str] = &[
    "SRSF_SR",
    "HNRNP",
    "U1_CORE",
    "U2_CORE",
    "SF3B_AXIS",
    "MINOR_U12",
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
    let mut entropy = vec![f32::NAN; n_cells];
    let mut dispersion = vec![f32::NAN; n_cells];

    let start = Instant::now();
    let ln_r = (n_reg as f32).ln();
    let mut zero_count = 0usize;

    for cell in 0..n_cells {
        let mut sum = 0.0f32;
        for &gene_id in &regulator_genes {
            sum += matrix.cp10k(gene_id as usize, cell);
        }
        if sum <= 0.0 {
            zero_count += 1;
            entropy[cell] = f32::NAN;
            dispersion[cell] = f32::NAN;
            continue;
        }

        let mut h = 0.0f32;
        let mut q = 0.0f32;
        for &gene_id in &regulator_genes {
            let y = matrix.cp10k(gene_id as usize, cell);
            let p = y / (sum + EPS);
            h -= p * (p + EPS).ln();
            q += p * p;
        }
        entropy[cell] = if ln_r > 0.0 { h / ln_r } else { f32::NAN };
        dispersion[cell] = (1.0 / (q + EPS)) / n_reg as f32;
    }

    if zero_count > 0 {
        warn!(
            zero_cells = zero_count,
            "cells with zero regulator expression"
        );
    }

    let med = median(&entropy);
    let mad_val = mad(&entropy, med) * 1.4826 + EPS;
    let mut z_entropy = vec![f32::NAN; n_cells];
    for cell in 0..n_cells {
        let value = entropy[cell];
        if value.is_finite() && med.is_finite() && mad_val.is_finite() {
            z_entropy[cell] = (value - med) / mad_val;
        } else {
            z_entropy[cell] = f32::NAN;
        }
    }

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
    ids_vec.sort();
    Ok(ids_vec)
}
