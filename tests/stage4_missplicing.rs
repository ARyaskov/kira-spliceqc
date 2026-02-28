use kira_spliceqc::model::geneset_activity::GenesetActivityMatrix;
use kira_spliceqc::pipeline::stage4_missplicing::compute;

fn make_activity(genesets: Vec<&str>, values: Vec<Vec<f32>>) -> GenesetActivityMatrix {
    let n_cells = values[0].len();
    let mut flat = Vec::new();
    for row in values {
        flat.extend_from_slice(&row);
    }
    GenesetActivityMatrix {
        genesets: genesets.into_iter().map(|s| s.to_string()).collect(),
        axes: vec!["AX".to_string(); flat.len() / n_cells],
        values: flat,
        n_cells,
    }
}

fn robust_z(values: &[f32]) -> Vec<f32> {
    let mut data: Vec<f32> = values.iter().copied().filter(|v| v.is_finite()).collect();
    data.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = if data.len() % 2 == 1 {
        data[data.len() / 2]
    } else {
        (data[data.len() / 2 - 1] + data[data.len() / 2]) * 0.5
    };
    let mut devs: Vec<f32> = values
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .map(|v| (v - med).abs())
        .collect();
    devs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mad = if devs.len() % 2 == 1 {
        devs[devs.len() / 2]
    } else {
        (devs[devs.len() / 2 - 1] + devs[devs.len() / 2]) * 0.5
    };
    let mad_scaled = mad * 1.4826 + 1e-6;
    values
        .iter()
        .map(|&v| {
            if v.is_finite() {
                (v - med) / mad_scaled
            } else {
                f32::NAN
            }
        })
        .collect()
}

#[test]
fn manual_z_scores_components() {
    let genesets = vec![
        "U1_CORE",
        "U2_CORE",
        "SF3B_AXIS",
        "SRSF_SR",
        "HNRNP",
        "MINOR_U12",
        "NMD_SURVEILLANCE",
    ];

    let values = vec![
        vec![-1.0, 1.0],
        vec![-1.0, 1.0],
        vec![-1.0, 1.0],
        vec![2.0, 0.0],
        vec![0.0, 2.0],
        vec![1.0, -1.0],
        vec![0.0, 2.0],
    ];

    let activity = make_activity(genesets, values);
    let metrics = compute(&activity).unwrap();

    let z_u1 = robust_z(&activity.values[0..2]);
    let z_u2 = robust_z(&activity.values[2..4]);
    let z_sf3b = robust_z(&activity.values[4..6]);
    let z_srsf = robust_z(&activity.values[6..8]);
    let z_hnrnp = robust_z(&activity.values[8..10]);
    let z_u12 = robust_z(&activity.values[10..12]);
    let z_nmd = robust_z(&activity.values[12..14]);

    let core_mean0 = (z_u1[0] + z_u2[0] + z_sf3b[0]) / 3.0;
    let b_core0 = (-core_mean0).max(0.0);
    let b_u120 = (z_u12[0] - (z_u1[0] + z_u2[0]) * 0.5).max(0.0);
    let b_nmd0 = z_nmd[0].max(0.0);
    let b_srhn0 = (z_srsf[0] - z_hnrnp[0]).abs();

    assert!((metrics.b_core[0] - b_core0).abs() < 1e-5);
    assert!((metrics.b_u12[0] - b_u120).abs() < 1e-5);
    assert!((metrics.b_nmd[0] - b_nmd0).abs() < 1e-5);
    assert!((metrics.b_srhn[0] - b_srhn0).abs() < 1e-5);

    let mb0 = metrics.burden[0];
    let mb_star0 = metrics.burden_star[0];
    let expected_mb0: f32 = 0.35 * b_core0 + 0.25 * b_u120 + 0.25 * b_nmd0 + 0.15 * b_srhn0;
    let expected_mb_star0 = 1.0_f32 - (-expected_mb0).exp();
    assert!((mb0 - expected_mb0).abs() < 1e-5);
    assert!((mb_star0 - expected_mb_star0).abs() < 1e-5);
}

#[test]
fn missing_geneset_nan_component() {
    let genesets = vec![
        "U1_CORE",
        "U2_CORE",
        "SF3B_AXIS",
        "SRSF_SR",
        "HNRNP",
        "MINOR_U12",
    ];
    let values = vec![
        vec![-1.0, 1.0],
        vec![-1.0, 1.0],
        vec![-1.0, 1.0],
        vec![2.0, 0.0],
        vec![0.0, 2.0],
        vec![1.0, -1.0],
    ];
    let activity = make_activity(genesets, values);
    let metrics = compute(&activity).unwrap();

    assert!(metrics.b_nmd[0].is_nan());
    assert!(metrics.b_core[0].is_finite());
    assert!(metrics.burden[0].is_nan());
}

#[test]
fn deterministic_outputs() {
    let genesets = vec![
        "U1_CORE",
        "U2_CORE",
        "SF3B_AXIS",
        "SRSF_SR",
        "HNRNP",
        "MINOR_U12",
        "NMD_SURVEILLANCE",
    ];
    let values = vec![
        vec![1.0, 2.0],
        vec![2.0, 3.0],
        vec![3.0, 4.0],
        vec![4.0, 5.0],
        vec![5.0, 6.0],
        vec![6.0, 7.0],
        vec![7.0, 8.0],
    ];
    let activity = make_activity(genesets, values);

    let first = compute(&activity).unwrap();
    let second = compute(&activity).unwrap();

    assert_eq!(first.b_core, second.b_core);
    assert_eq!(first.b_u12, second.b_u12);
    assert_eq!(first.b_nmd, second.b_nmd);
    assert_eq!(first.b_srhn, second.b_srhn);
    assert_eq!(first.burden, second.burden);
    assert_eq!(first.burden_star, second.burden_star);
}
