use rayon::prelude::*;
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

    let rows: Vec<(f32, SpliceIntegrityClass, f32, f32, f32, f32)> = (0..n_cells)
        .into_par_iter()
        .map(|cell| {
            let mb_star = missplicing.burden_star[cell];
            let i_val = imbalance.imbalance[cell];
            let z_ent = isoform.z_entropy[cell];
            let ent = isoform.entropy[cell];

            let p1 = mb_star;
            let p2 = clamp01((i_val - 0.8) / 1.2);
            let p3 = clamp01((z_ent - 1.5) / 2.0);
            let p4 = clamp01((ent - 0.85) / 0.15);

            if [p1, p2, p3, p4].iter().all(|v| v.is_finite()) {
                let raw = 1.0 - (0.35 * p1 + 0.25 * p2 + 0.25 * p3 + 0.15 * p4);
                let score = clamp01(raw);
                (score, classify(score), p1, p2, p3, p4)
            } else {
                (f32::NAN, SpliceIntegrityClass::Broken, p1, p2, p3, p4)
            }
        })
        .collect();

    let mut sis = Vec::with_capacity(n_cells);
    let mut class = Vec::with_capacity(n_cells);
    let mut p_missplicing = Vec::with_capacity(n_cells);
    let mut p_imbalance = Vec::with_capacity(n_cells);
    let mut p_entropy_z = Vec::with_capacity(n_cells);
    let mut p_entropy_abs = Vec::with_capacity(n_cells);

    for (s, c, p1, p2, p3, p4) in rows {
        sis.push(s);
        class.push(c);
        p_missplicing.push(p1);
        p_imbalance.push(p2);
        p_entropy_z.push(p3);
        p_entropy_abs.push(p4);
    }

    let mut min = f32::INFINITY;
    let mut max = f32::NEG_INFINITY;
    for &v in &sis {
        if v.is_finite() {
            if v < min {
                min = v;
            }
            if v > max {
                max = v;
            }
        }
    }
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

#[inline]
fn clamp01(value: f32) -> f32 {
    if !value.is_finite() {
        f32::NAN
    } else {
        value.clamp(0.0, 1.0)
    }
}

#[inline]
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
