use kira_spliceqc::model::geneset_activity::GenesetActivityMatrix;
use kira_spliceqc::pipeline::stage8_coupling::compute;

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
fn synthetic_coupling_stress() {
    let genesets = vec!["TRANSCRIPTION_COUPLING", "U1_CORE", "U2_CORE", "SF3B_AXIS"];
    let values = vec![
        vec![2.0, 0.0],
        vec![0.0, 2.0],
        vec![0.0, 2.0],
        vec![0.0, 2.0],
    ];
    let activity = make_activity(genesets, values);
    let metrics = compute(&activity).unwrap();

    let z_tx = robust_z(&activity.values[0..2]);
    let z_u1 = robust_z(&activity.values[2..4]);
    let z_u2 = robust_z(&activity.values[4..6]);
    let z_sf3b = robust_z(&activity.values[6..8]);

    let expected0 = z_tx[0] - (z_u1[0] + z_u2[0] + z_sf3b[0]) / 3.0;
    let expected1 = z_tx[1] - (z_u1[1] + z_u2[1] + z_sf3b[1]) / 3.0;

    assert!((metrics.coupling_stress[0] - expected0).abs() < 1e-5);
    assert!((metrics.coupling_stress[1] - expected1).abs() < 1e-5);
}

#[test]
fn missing_coupling_geneset_errors() {
    let genesets = vec!["U1_CORE", "U2_CORE", "SF3B_AXIS"];
    let values = vec![vec![0.0, 1.0], vec![0.0, 1.0], vec![0.0, 1.0]];
    let activity = make_activity(genesets, values);
    let err = compute(&activity).unwrap_err();
    assert!(format!("{err}").contains("coupling"));
}

#[test]
fn deterministic_outputs() {
    let genesets = vec!["TRANSCRIPTION_COUPLING", "U1_CORE", "U2_CORE", "SF3B_AXIS"];
    let values = vec![
        vec![1.0, 2.0],
        vec![2.0, 3.0],
        vec![3.0, 4.0],
        vec![4.0, 5.0],
    ];
    let activity = make_activity(genesets, values);
    let first = compute(&activity).unwrap();
    let second = compute(&activity).unwrap();

    assert_eq!(first.coupling_stress, second.coupling_stress);
}
