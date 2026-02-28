use kira_spliceqc::model::timecourse::SplicingTrajectoryClass;
use kira_spliceqc::pipeline::stage14_timecourse::compute;

#[test]
fn synthetic_adaptive_trajectory() {
    let timepoints = vec![0, 1, 2];
    let sis = vec![0.5, 0.6, 0.7];
    let entropy = vec![2.0, 1.5, 1.2];
    let imbalance = vec![1.5, 1.2, 1.0];

    let metrics = compute(&timepoints, &sis, &entropy, &imbalance).unwrap();
    assert_eq!(metrics.trajectory, SplicingTrajectoryClass::Adaptive);
}

#[test]
fn synthetic_degenerative_trajectory() {
    let timepoints = vec![0, 1, 2];
    let sis = vec![0.7, 0.6, 0.4];
    let entropy = vec![1.0, 1.5, 2.0];
    let imbalance = vec![1.0, 1.8, 2.2];

    let metrics = compute(&timepoints, &sis, &entropy, &imbalance).unwrap();
    assert_eq!(metrics.trajectory, SplicingTrajectoryClass::Degenerative);
}

#[test]
fn oscillatory_pattern() {
    let timepoints = vec![0, 1, 2];
    let sis = vec![0.5, 0.6, 0.55];
    let entropy = vec![1.0, 1.1, 1.0];
    let imbalance = vec![1.0, 0.9, 1.1];

    let metrics = compute(&timepoints, &sis, &entropy, &imbalance).unwrap();
    assert_eq!(metrics.trajectory, SplicingTrajectoryClass::Oscillatory);
}

#[test]
fn too_few_timepoints_inconclusive() {
    let timepoints = vec![0, 1];
    let sis = vec![0.5, 0.6];
    let entropy = vec![1.0, 0.9];
    let imbalance = vec![1.0, 0.8];

    let metrics = compute(&timepoints, &sis, &entropy, &imbalance).unwrap();
    assert_eq!(metrics.trajectory, SplicingTrajectoryClass::Inconclusive);
}

#[test]
fn deterministic_outputs() {
    let timepoints = vec![0, 1, 2, 3];
    let sis = vec![0.5, 0.6, 0.7, 0.8];
    let entropy = vec![1.5, 1.4, 1.3, 1.2];
    let imbalance = vec![1.0, 0.9, 0.8, 0.7];

    let first = compute(&timepoints, &sis, &entropy, &imbalance).unwrap();
    let second = compute(&timepoints, &sis, &entropy, &imbalance).unwrap();

    assert_eq!(first.trajectory, second.trajectory);
    assert_eq!(first.delta_sis, second.delta_sis);
    assert_eq!(first.delta_entropy, second.delta_entropy);
    assert_eq!(first.delta_imbalance, second.delta_imbalance);
}
