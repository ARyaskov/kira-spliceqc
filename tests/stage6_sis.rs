use kira_spliceqc::model::imbalance::SpliceosomeImbalanceMetrics;
use kira_spliceqc::model::isoform_dispersion::IsoformDispersionMetrics;
use kira_spliceqc::model::missplicing::MissplicingMetrics;
use kira_spliceqc::model::sis::SpliceIntegrityClass;
use kira_spliceqc::pipeline::stage6_sis::run_stage6;

fn make_isoform(entropy: Vec<f32>, z_entropy: Vec<f32>) -> IsoformDispersionMetrics {
    IsoformDispersionMetrics {
        entropy,
        dispersion: Vec::new(),
        z_entropy,
    }
}

fn make_missplicing(burden_star: Vec<f32>) -> MissplicingMetrics {
    MissplicingMetrics {
        b_core: Vec::new(),
        b_u12: Vec::new(),
        b_nmd: Vec::new(),
        b_srhn: Vec::new(),
        burden: Vec::new(),
        burden_star,
    }
}

fn make_imbalance(imbalance: Vec<f32>) -> SpliceosomeImbalanceMetrics {
    let n = imbalance.len();
    SpliceosomeImbalanceMetrics {
        z_u1: vec![0.0; n],
        z_u2: vec![0.0; n],
        z_sf3b: vec![0.0; n],
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

#[test]
fn synthetic_penalties_exact_sis() {
    let isoform = make_isoform(vec![0.90], vec![2.5]);
    let missplicing = make_missplicing(vec![0.2]);
    let imbalance = make_imbalance(vec![1.4]);

    let metrics = run_stage6(&isoform, &missplicing, &imbalance).unwrap();

    let p1 = 0.2_f32;
    let p2 = ((1.4_f32 - 0.8) / 1.2).clamp(0.0, 1.0);
    let p3 = ((2.5_f32 - 1.5) / 2.0).clamp(0.0, 1.0);
    let p4 = ((0.90_f32 - 0.85) / 0.15).clamp(0.0, 1.0);
    let expected = 1.0 - (0.35 * p1 + 0.25 * p2 + 0.25 * p3 + 0.15 * p4);
    let expected = expected.clamp(0.0, 1.0);

    assert!((metrics.sis[0] - expected).abs() < 1e-6);
}

#[test]
fn class_boundaries() {
    let isoform = make_isoform(vec![0.85; 4], vec![0.0; 4]);
    let missplicing = make_missplicing(vec![0.0; 4]);
    let imbalance = make_imbalance(vec![0.0; 4]);
    let metrics = run_stage6(&isoform, &missplicing, &imbalance).unwrap();

    assert_eq!(metrics.class[0], SpliceIntegrityClass::Intact);

    // Construct penalties within [0,1] to hit boundary classes
    // Cell1: p1=0.6 -> SIS=0.79 (Stressed)
    // Cell2: p1=1.0, p2=0.2 -> SIS=0.60 (Stressed boundary)
    // Cell3: p1=1.0, p2=1.0, p3=0.4 -> SIS=0.30 (Broken)
    let missplicing = make_missplicing(vec![0.0, 0.6, 1.0, 1.0]);
    let isoform = make_isoform(vec![0.85; 4], vec![1.5, 1.5, 1.5, 2.3]);
    let imbalance = make_imbalance(vec![0.0, 0.0, 1.04, 2.0]);
    let metrics = run_stage6(&isoform, &missplicing, &imbalance).unwrap();

    assert_eq!(metrics.class[0], SpliceIntegrityClass::Intact);
    assert_eq!(metrics.class[1], SpliceIntegrityClass::Stressed);
    assert_eq!(metrics.class[2], SpliceIntegrityClass::Stressed);
    assert_eq!(metrics.class[3], SpliceIntegrityClass::Broken);
}

#[test]
fn nan_propagation() {
    let isoform = make_isoform(vec![f32::NAN], vec![0.0]);
    let missplicing = make_missplicing(vec![0.0]);
    let imbalance = make_imbalance(vec![0.0]);
    let metrics = run_stage6(&isoform, &missplicing, &imbalance).unwrap();

    assert!(metrics.sis[0].is_nan());
    assert_eq!(metrics.class[0], SpliceIntegrityClass::Broken);
}

#[test]
fn deterministic_outputs() {
    let isoform = make_isoform(vec![0.9, 0.8], vec![1.0, 2.0]);
    let missplicing = make_missplicing(vec![0.1, 0.2]);
    let imbalance = make_imbalance(vec![1.0, 1.5]);

    let first = run_stage6(&isoform, &missplicing, &imbalance).unwrap();
    let second = run_stage6(&isoform, &missplicing, &imbalance).unwrap();

    assert_eq!(first.sis, second.sis);
    assert_eq!(first.class, second.class);
    assert_eq!(first.p_missplicing, second.p_missplicing);
    assert_eq!(first.p_imbalance, second.p_imbalance);
    assert_eq!(first.p_entropy_z, second.p_entropy_z);
    assert_eq!(first.p_entropy_abs, second.p_entropy_abs);
}
