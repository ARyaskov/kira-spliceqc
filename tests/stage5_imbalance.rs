use kira_spliceqc::model::geneset_activity::GenesetActivityMatrix;
use kira_spliceqc::pipeline::stage5_imbalance::compute;

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
fn synthetic_z_axes_and_magnitude() {
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
        vec![1.0, -1.0],
        vec![0.0, 0.0],
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

    assert!((metrics.axis_sr_hnrnp[0] - (z_srsf[0] - z_hnrnp[0])).abs() < 1e-5);
    assert!((metrics.axis_u2_u1[0] - (z_u2[0] - z_u1[0])).abs() < 1e-5);
    assert!((metrics.axis_u12_major[0] - (z_u12[0] - (z_u1[0] + z_u2[0]) * 0.5)).abs() < 1e-5);
    assert!((metrics.axis_nmd[0] - z_nmd[0]).abs() < 1e-5);

    let mut vals = [z_u1[0], z_u2[0], z_sf3b[0], z_srsf[0], z_hnrnp[0], z_u12[0]];
    for v in &mut vals {
        if *v < -6.0 {
            *v = -6.0;
        } else if *v > 6.0 {
            *v = 6.0;
        }
    }
    let mean_sq = vals.iter().map(|v| v * v).sum::<f32>() / vals.len() as f32;
    let expected_imb = mean_sq.sqrt();
    assert!((metrics.imbalance[0] - expected_imb).abs() < 1e-5);
}

#[test]
fn missing_geneset_nan_propagates() {
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
        vec![1.0, -1.0],
        vec![0.0, 0.0],
        vec![2.0, 0.0],
        vec![0.0, 2.0],
        vec![1.0, -1.0],
    ];
    let activity = make_activity(genesets, values);
    let metrics = compute(&activity).unwrap();

    assert!(metrics.axis_nmd[0].is_nan());
    assert!(metrics.imbalance[0].is_finite());
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

    assert_eq!(first.z_u1, second.z_u1);
    assert_eq!(first.z_u2, second.z_u2);
    assert_eq!(first.z_sf3b, second.z_sf3b);
    assert_eq!(first.z_srsf, second.z_srsf);
    assert_eq!(first.z_hnrnp, second.z_hnrnp);
    assert_eq!(first.z_u12, second.z_u12);
    assert_eq!(first.z_nmd, second.z_nmd);
    assert_eq!(first.axis_sr_hnrnp, second.axis_sr_hnrnp);
    assert_eq!(first.axis_u2_u1, second.axis_u2_u1);
    assert_eq!(first.axis_u12_major, second.axis_u12_major);
    assert_eq!(first.axis_nmd, second.axis_nmd);
    assert_eq!(first.imbalance, second.imbalance);
}
