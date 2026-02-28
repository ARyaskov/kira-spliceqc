use kira_spliceqc::model::cryptic_risk::CrypticSplicingRiskMetrics;
use kira_spliceqc::model::imbalance::SpliceosomeImbalanceMetrics;
use kira_spliceqc::model::isoform_dispersion::IsoformDispersionMetrics;
use kira_spliceqc::pipeline::stage12_cryptic_risk::compute;

fn make_isoform(z_entropy: Vec<f32>) -> IsoformDispersionMetrics {
    IsoformDispersionMetrics {
        entropy: vec![0.0; z_entropy.len()],
        dispersion: vec![0.0; z_entropy.len()],
        z_entropy,
    }
}

fn make_imbalance(axis_sr_hnrnp: Vec<f32>, z_nmd: Vec<f32>) -> SpliceosomeImbalanceMetrics {
    let n = axis_sr_hnrnp.len();
    SpliceosomeImbalanceMetrics {
        z_u1: vec![0.0; n],
        z_u2: vec![0.0; n],
        z_sf3b: vec![0.0; n],
        z_srsf: vec![0.0; n],
        z_hnrnp: vec![0.0; n],
        z_u12: vec![0.0; n],
        z_nmd,
        axis_sr_hnrnp,
        axis_u2_u1: vec![0.0; n],
        axis_u12_major: vec![0.0; n],
        axis_nmd: vec![0.0; n],
        imbalance: vec![0.0; n],
    }
}

fn risk(metrics: &CrypticSplicingRiskMetrics, cell: usize) -> f32 {
    metrics.cryptic_risk[cell]
}

#[test]
fn synthetic_signal_combinations() {
    let isoform = make_isoform(vec![0.0, 6.0]);
    let imbalance = make_imbalance(vec![6.0, 6.0], vec![0.0, 6.0]);
    let metrics = compute(&isoform, &imbalance).unwrap();

    let single = risk(&metrics, 0);
    let double = risk(&metrics, 1);

    assert!(single > 0.30 && single < 0.45);
    assert!(double > 0.55);
}

#[test]
fn nan_propagation() {
    let isoform = make_isoform(vec![1.0, 1.0]);
    let imbalance = make_imbalance(vec![f32::NAN, 1.0], vec![1.0, 1.0]);
    let metrics = compute(&isoform, &imbalance).unwrap();

    assert!(metrics.cryptic_risk[0].is_nan());
    assert!(metrics.x_sr_hnrnp[0].is_nan());
    assert!(metrics.x_entropy[0].is_nan());
    assert!(metrics.x_nmd[0].is_nan());
    assert!(metrics.cryptic_risk[1].is_finite());
}

#[test]
fn boundary_behavior() {
    let isoform = make_isoform(vec![0.0, 6.0]);
    let imbalance = make_imbalance(vec![0.0, 6.0], vec![0.0, 6.0]);
    let metrics = compute(&isoform, &imbalance).unwrap();

    let low = risk(&metrics, 0);
    let high = risk(&metrics, 1);

    assert!((low - 0.182425).abs() < 1e-3);
    assert!((high - 0.817574).abs() < 1e-3);
}

#[test]
fn deterministic_outputs() {
    let isoform = make_isoform(vec![1.0, 2.0, 3.0]);
    let imbalance = make_imbalance(vec![2.0, 3.0, 4.0], vec![3.0, 4.0, 5.0]);

    let first = compute(&isoform, &imbalance).unwrap();
    let second = compute(&isoform, &imbalance).unwrap();

    assert_eq!(first.cryptic_risk, second.cryptic_risk);
    assert_eq!(first.x_sr_hnrnp, second.x_sr_hnrnp);
    assert_eq!(first.x_entropy, second.x_entropy);
    assert_eq!(first.x_nmd, second.x_nmd);
}
