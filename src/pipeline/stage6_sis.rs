use tracing::{debug, info};

use crate::input::error::InputError;
use crate::model::imbalance::SpliceosomeImbalanceMetrics;
use crate::model::isoform_dispersion::IsoformDispersionMetrics;
use crate::model::missplicing::MissplicingMetrics;
use crate::model::sis::{SpliceIntegrityClass, SpliceIntegrityMetrics};
use crate::stats::robust::median;

pub fn run_stage6(
    isoform: &IsoformDispersionMetrics,
    missplicing: &MissplicingMetrics,
    imbalance: &SpliceosomeImbalanceMetrics,
) -> Result<SpliceIntegrityMetrics, InputError> {
    info!("SIS computation started");
    validate_lengths(isoform, missplicing, imbalance)?;

    let n_cells = isoform.entropy.len();
    let mut sis = vec![f32::NAN; n_cells];
    let mut class = vec![SpliceIntegrityClass::Broken; n_cells];
    let mut p_missplicing = vec![f32::NAN; n_cells];
    let mut p_imbalance = vec![f32::NAN; n_cells];
    let mut p_entropy_z = vec![f32::NAN; n_cells];
    let mut p_entropy_abs = vec![f32::NAN; n_cells];

    for cell in 0..n_cells {
        let mb_star = missplicing.burden_star[cell];
        let i_val = imbalance.imbalance[cell];
        let z_ent = isoform.z_entropy[cell];
        let ent = isoform.entropy[cell];

        let p1 = mb_star;
        let p2 = clamp01((i_val - 0.8) / 1.2);
        let p3 = clamp01((z_ent - 1.5) / 2.0);
        let p4 = clamp01((ent - 0.85) / 0.15);

        p_missplicing[cell] = p1;
        p_imbalance[cell] = p2;
        p_entropy_z[cell] = p3;
        p_entropy_abs[cell] = p4;

        if [p1, p2, p3, p4].iter().all(|v| v.is_finite()) {
            let raw = 1.0 - (0.35 * p1 + 0.25 * p2 + 0.25 * p3 + 0.15 * p4);
            let score = clamp01(raw);
            sis[cell] = score;
            class[cell] = classify(score);
        } else {
            sis[cell] = f32::NAN;
            class[cell] = SpliceIntegrityClass::Broken;
        }
    }

    let min = sis
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .fold(f32::INFINITY, f32::min);
    let max = sis
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .fold(f32::NEG_INFINITY, f32::max);
    let med = median(&sis);
    debug!(min = min, median = med, max = max, "SIS summary");

    Ok(SpliceIntegrityMetrics {
        sis,
        class,
        p_missplicing,
        p_imbalance,
        p_entropy_z,
        p_entropy_abs,
    })
}

fn validate_lengths(
    isoform: &IsoformDispersionMetrics,
    missplicing: &MissplicingMetrics,
    imbalance: &SpliceosomeImbalanceMetrics,
) -> Result<(), InputError> {
    let n = isoform.entropy.len();
    if isoform.z_entropy.len() != n {
        return Err(InputError::LengthMismatch(
            "isoform entropy/z_entropy length mismatch".to_string(),
        ));
    }
    if missplicing.burden_star.len() != n {
        return Err(InputError::LengthMismatch(
            "missplicing length mismatch".to_string(),
        ));
    }
    if imbalance.imbalance.len() != n {
        return Err(InputError::LengthMismatch(
            "imbalance length mismatch".to_string(),
        ));
    }
    Ok(())
}

fn clamp01(value: f32) -> f32 {
    if !value.is_finite() {
        f32::NAN
    } else if value < 0.0 {
        0.0
    } else if value > 1.0 {
        1.0
    } else {
        value
    }
}

fn classify(score: f32) -> SpliceIntegrityClass {
    if score >= 0.80 {
        SpliceIntegrityClass::Intact
    } else if score >= 0.60 {
        SpliceIntegrityClass::Stressed
    } else if score >= 0.40 {
        SpliceIntegrityClass::Impaired
    } else {
        SpliceIntegrityClass::Broken
    }
}
