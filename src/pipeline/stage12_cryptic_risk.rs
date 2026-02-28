use tracing::{debug, info};

use crate::input::error::InputError;
use crate::model::cryptic_risk::CrypticSplicingRiskMetrics;
use crate::model::imbalance::SpliceosomeImbalanceMetrics;
use crate::model::isoform_dispersion::IsoformDispersionMetrics;
use crate::stats::robust::median;

const SAT: f32 = 6.0;

pub fn run_stage12(
    isoform: &IsoformDispersionMetrics,
    imbalance: &SpliceosomeImbalanceMetrics,
) -> Result<CrypticSplicingRiskMetrics, InputError> {
    compute(isoform, imbalance)
}

pub fn compute(
    isoform: &IsoformDispersionMetrics,
    imbalance: &SpliceosomeImbalanceMetrics,
) -> Result<CrypticSplicingRiskMetrics, InputError> {
    info!("cryptic splicing risk computation started");
    validate_lengths(isoform, imbalance)?;

    let n_cells = isoform.z_entropy.len();
    if n_cells == 0 {
        return Err(InputError::MissingCrypticRiskSignals);
    }

    let mut cryptic_risk = vec![f32::NAN; n_cells];
    let mut x_sr_hnrnp = vec![f32::NAN; n_cells];
    let mut x_entropy = vec![f32::NAN; n_cells];
    let mut x_nmd = vec![f32::NAN; n_cells];

    for cell in 0..n_cells {
        let axis = imbalance.axis_sr_hnrnp[cell];
        let z_ent = isoform.z_entropy[cell];
        let z_nmd = imbalance.z_nmd[cell];

        if axis.is_finite() && z_ent.is_finite() && z_nmd.is_finite() {
            let x1 = clamp01(axis.abs());
            let x2 = clamp01(z_ent);
            let x3 = clamp01(z_nmd);

            x_sr_hnrnp[cell] = x1;
            x_entropy[cell] = x2;
            x_nmd[cell] = x3;

            let raw = x1 + x2 + x3;
            cryptic_risk[cell] = sigmoid(raw - 1.5);
        } else {
            cryptic_risk[cell] = f32::NAN;
            x_sr_hnrnp[cell] = f32::NAN;
            x_entropy[cell] = f32::NAN;
            x_nmd[cell] = f32::NAN;
        }
    }

    let (min, max) = finite_min_max(&cryptic_risk);
    let med = median(&cryptic_risk);
    debug!(min = min, median = med, max = max, "cryptic risk summary");

    Ok(CrypticSplicingRiskMetrics {
        cryptic_risk,
        x_sr_hnrnp,
        x_entropy,
        x_nmd,
    })
}

fn validate_lengths(
    isoform: &IsoformDispersionMetrics,
    imbalance: &SpliceosomeImbalanceMetrics,
) -> Result<(), InputError> {
    let n = isoform.z_entropy.len();
    if isoform.entropy.len() != n {
        return Err(InputError::LengthMismatch(
            "isoform entropy/z_entropy length mismatch".to_string(),
        ));
    }
    if imbalance.axis_sr_hnrnp.len() != n {
        return Err(InputError::LengthMismatch(
            "imbalance axis_sr_hnrnp length mismatch".to_string(),
        ));
    }
    if imbalance.z_nmd.len() != n {
        return Err(InputError::LengthMismatch(
            "imbalance z_nmd length mismatch".to_string(),
        ));
    }
    Ok(())
}

fn clamp01(value: f32) -> f32 {
    if !value.is_finite() {
        f32::NAN
    } else if value <= 0.0 {
        0.0
    } else if value >= SAT {
        1.0
    } else {
        value / SAT
    }
}

fn sigmoid(value: f32) -> f32 {
    1.0 / (1.0 + (-value).exp())
}

fn finite_min_max(values: &[f32]) -> (f32, f32) {
    let mut min = f32::INFINITY;
    let mut max = f32::NEG_INFINITY;
    let mut found = false;
    for &v in values {
        if v.is_finite() {
            found = true;
            if v < min {
                min = v;
            }
            if v > max {
                max = v;
            }
        }
    }
    if found {
        (min, max)
    } else {
        (f32::NAN, f32::NAN)
    }
}
