use kira_spliceqc::model::collapse::SpliceosomeCollapseStatus;
use kira_spliceqc::model::imbalance::SpliceosomeImbalanceMetrics;
use kira_spliceqc::model::sis::{SpliceIntegrityClass, SpliceIntegrityMetrics};
use kira_spliceqc::pipeline::stage13_collapse::compute;

fn make_imbalance(
    z_u1: Vec<f32>,
    z_u2: Vec<f32>,
    z_sf3b: Vec<f32>,
    imbalance: Vec<f32>,
) -> SpliceosomeImbalanceMetrics {
    let n = z_u1.len();
    SpliceosomeImbalanceMetrics {
        z_u1,
        z_u2,
        z_sf3b,
        z_srsf: vec![0.0; n],
        z_hnrnp: vec![0.0; n],
        z_u12: vec![0.0; n],
        z_nmd: vec![0.0; n],
        axis_sr_hnrnp: vec![0.0; n],
        axis_u2_u1: vec![0.0; n],
        axis_u12_major: vec![0.0; n],
        axis_nmd: vec![0.0; n],
        imbalance,
    }
}

fn make_sis(values: Vec<f32>) -> SpliceIntegrityMetrics {
    let n = values.len();
    SpliceIntegrityMetrics {
        sis: values,
        class: vec![SpliceIntegrityClass::Broken; n],
        p_missplicing: vec![0.0; n],
        p_imbalance: vec![0.0; n],
        p_entropy_z: vec![0.0; n],
        p_entropy_abs: vec![0.0; n],
    }
}

#[test]
fn collapsed_cell_detected() {
    let imbalance = make_imbalance(vec![-2.0], vec![-2.1], vec![-1.6], vec![2.0]);
    let sis = make_sis(vec![0.2]);
    let metrics = compute(&imbalance, &sis).unwrap();

    assert_eq!(
        metrics.collapse_status[0],
        SpliceosomeCollapseStatus::Collapse
    );
    assert!(metrics.core_suppression[0]);
    assert!(metrics.high_imbalance[0]);
    assert!(metrics.low_sis[0]);
}

#[test]
fn partial_failure_is_no_collapse() {
    let imbalance = make_imbalance(vec![-2.0], vec![-2.1], vec![-1.6], vec![1.0]);
    let sis = make_sis(vec![0.2]);
    let metrics = compute(&imbalance, &sis).unwrap();

    assert_eq!(
        metrics.collapse_status[0],
        SpliceosomeCollapseStatus::NoCollapse
    );
    assert!(metrics.core_suppression[0]);
    assert!(!metrics.high_imbalance[0]);
    assert!(metrics.low_sis[0]);
}

#[test]
fn nan_core_z_is_inconclusive() {
    let imbalance = make_imbalance(vec![f32::NAN], vec![-2.0], vec![-2.0], vec![2.0]);
    let sis = make_sis(vec![0.2]);
    let metrics = compute(&imbalance, &sis).unwrap();

    assert_eq!(
        metrics.collapse_status[0],
        SpliceosomeCollapseStatus::Inconclusive
    );
}

#[test]
fn deterministic_outputs() {
    let imbalance = make_imbalance(
        vec![-2.0, 0.0],
        vec![-2.0, 0.0],
        vec![-2.0, 0.0],
        vec![2.0, 0.5],
    );
    let sis = make_sis(vec![0.2, 0.9]);

    let first = compute(&imbalance, &sis).unwrap();
    let second = compute(&imbalance, &sis).unwrap();

    assert_eq!(first.collapse_status, second.collapse_status);
    assert_eq!(first.core_suppression, second.core_suppression);
    assert_eq!(first.high_imbalance, second.high_imbalance);
    assert_eq!(first.low_sis, second.low_sis);
}
