use rayon::prelude::*;
use tracing::{debug, info};

use crate::input::error::InputError;
use crate::model::collapse::{SpliceosomeCollapseMetrics, SpliceosomeCollapseStatus};
use crate::model::imbalance::SpliceosomeImbalanceMetrics;
use crate::model::sis::SpliceIntegrityMetrics;

pub fn run_stage13(
    imbalance: &SpliceosomeImbalanceMetrics,
    sis: &SpliceIntegrityMetrics,
) -> Result<SpliceosomeCollapseMetrics, InputError> {
    compute(imbalance, sis)
}

pub fn compute(
    imbalance: &SpliceosomeImbalanceMetrics,
    sis: &SpliceIntegrityMetrics,
) -> Result<SpliceosomeCollapseMetrics, InputError> {
    validate_lengths(imbalance, sis)?;

    let n_cells = imbalance.z_u1.len();
    if n_cells == 0 {
        return Err(InputError::MissingCollapseSignals);
    }

    let rows: Vec<(SpliceosomeCollapseStatus, bool, bool, bool)> = (0..n_cells)
        .into_par_iter()
        .map(|cell| {
            let z_u1 = imbalance.z_u1[cell];
            let z_u2 = imbalance.z_u2[cell];
            let z_sf3b = imbalance.z_sf3b[cell];
            let i_val = imbalance.imbalance[cell];
            let sis_val = sis.sis[cell];

            let core_missing = !z_u1.is_finite() || !z_u2.is_finite() || !z_sf3b.is_finite();
            let high_imbalance = i_val.is_finite() && i_val > 1.8;
            let low_sis = sis_val.is_finite() && sis_val < 0.4;

            if core_missing {
                return (
                    SpliceosomeCollapseStatus::Inconclusive,
                    false,
                    high_imbalance,
                    low_sis,
                );
            }

            let core_suppression = z_u1 < -1.5 && z_u2 < -1.5 && z_sf3b < -1.5;
            let status = if core_suppression && high_imbalance && low_sis {
                SpliceosomeCollapseStatus::Collapse
            } else {
                SpliceosomeCollapseStatus::NoCollapse
            };
            (status, core_suppression, high_imbalance, low_sis)
        })
        .collect();

    let mut collapse_status = Vec::with_capacity(n_cells);
    let mut core_suppression = Vec::with_capacity(n_cells);
    let mut high_imbalance = Vec::with_capacity(n_cells);
    let mut low_sis = Vec::with_capacity(n_cells);

    let mut count_collapse = 0usize;
    let mut count_no = 0usize;
    let mut count_inconclusive = 0usize;
    for (s, c, h, l) in rows {
        match s {
            SpliceosomeCollapseStatus::Collapse => count_collapse += 1,
            SpliceosomeCollapseStatus::NoCollapse => count_no += 1,
            SpliceosomeCollapseStatus::Inconclusive => count_inconclusive += 1,
        }
        collapse_status.push(s);
        core_suppression.push(c);
        high_imbalance.push(h);
        low_sis.push(l);
    }

    info!("collapse detection completed");
    debug!(
        collapse = count_collapse,
        no_collapse = count_no,
        inconclusive = count_inconclusive,
        "collapse status counts"
    );

    Ok(SpliceosomeCollapseMetrics {
        collapse_status,
        core_suppression,
        high_imbalance,
        low_sis,
    })
}

fn validate_lengths(
    imbalance: &SpliceosomeImbalanceMetrics,
    sis: &SpliceIntegrityMetrics,
) -> Result<(), InputError> {
    let n = imbalance.z_u1.len();
    if imbalance.z_u2.len() != n || imbalance.z_sf3b.len() != n {
        return Err(InputError::LengthMismatch(
            "imbalance core z-score length mismatch".to_string(),
        ));
    }
    if imbalance.imbalance.len() != n {
        return Err(InputError::LengthMismatch(
            "imbalance integrated length mismatch".to_string(),
        ));
    }
    if sis.sis.len() != n {
        return Err(InputError::LengthMismatch(
            "SIS length mismatch".to_string(),
        ));
    }
    Ok(())
}
