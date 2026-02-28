use tracing::{debug, info};

use crate::input::error::InputError;
use crate::model::timecourse::{SplicingTrajectoryClass, TimecourseSplicingMetrics};

pub fn run_stage14(
    timepoints: Option<&[i64]>,
    sis_median: &[f32],
    entropy_median: &[f32],
    imbalance_median: &[f32],
) -> Result<Option<TimecourseSplicingMetrics>, InputError> {
    let Some(timepoints) = timepoints else {
        return Ok(None);
    };
    let metrics = compute(timepoints, sis_median, entropy_median, imbalance_median)?;
    Ok(Some(metrics))
}

pub fn compute(
    timepoints: &[i64],
    sis_median: &[f32],
    entropy_median: &[f32],
    imbalance_median: &[f32],
) -> Result<TimecourseSplicingMetrics, InputError> {
    let n = timepoints.len();
    if n == 0 {
        return Err(InputError::MissingTimecourseSignals);
    }
    if sis_median.len() != n || entropy_median.len() != n || imbalance_median.len() != n {
        return Err(InputError::LengthMismatch(
            "timecourse median length mismatch".to_string(),
        ));
    }

    let delta_sis = deltas(sis_median);
    let delta_entropy = deltas(entropy_median);
    let delta_imbalance = deltas(imbalance_median);

    let mut trajectory = SplicingTrajectoryClass::Inconclusive;
    let has_nan = sis_median
        .iter()
        .chain(entropy_median.iter())
        .chain(imbalance_median.iter())
        .any(|v| !v.is_finite())
        || delta_sis
            .iter()
            .chain(delta_entropy.iter())
            .chain(delta_imbalance.iter())
            .any(|v| !v.is_finite());

    if n >= 3 && !has_nan {
        let adaptive = delta_sis.iter().all(|d| *d > 0.0)
            && delta_entropy.iter().all(|d| *d <= 0.0)
            && delta_imbalance.iter().all(|d| *d <= 0.0);
        let degenerative = delta_sis.iter().all(|d| *d < 0.0)
            && delta_entropy.iter().all(|d| *d >= 0.0)
            && delta_imbalance.iter().all(|d| *d >= 0.0);

        trajectory = if adaptive {
            SplicingTrajectoryClass::Adaptive
        } else if degenerative {
            SplicingTrajectoryClass::Degenerative
        } else {
            SplicingTrajectoryClass::Oscillatory
        };
    }

    info!("time-course coherence analysis completed");
    debug!(
        trajectory = ?trajectory,
        delta_sis = ?delta_sis,
        delta_entropy = ?delta_entropy,
        delta_imbalance = ?delta_imbalance,
        "time-course deltas"
    );

    Ok(TimecourseSplicingMetrics {
        trajectory,
        delta_sis,
        delta_entropy,
        delta_imbalance,
    })
}

fn deltas(values: &[f32]) -> Vec<f32> {
    if values.len() < 2 {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(values.len() - 1);
    for idx in 0..(values.len() - 1) {
        out.push(values[idx + 1] - values[idx]);
    }
    out
}
