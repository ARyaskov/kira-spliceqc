use kira_spliceqc::model::geneset_activity::GenesetActivityMatrix;
use kira_spliceqc::pipeline::stage10_assembly_phase::compute;

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
fn synthetic_phase_imbalance() {
    let genesets = vec![
        "SPLICE_EA_PHASE",
        "SPLICE_B_PHASE",
        "SPLICE_CATALYTIC_PHASE",
    ];
    let values = vec![vec![2.0, 0.0], vec![0.0, 2.0], vec![1.0, 1.0]];
    let activity = make_activity(genesets, values);
    let metrics = compute(&activity).unwrap();

    let z_ea = robust_z(&activity.values[0..2]);
    let z_b = robust_z(&activity.values[2..4]);
    let z_cat = robust_z(&activity.values[4..6]);

    let expected_ea0 = z_ea[0] - (z_b[0] + z_cat[0]) * 0.5;
    let expected_b0 = z_b[0] - (z_ea[0] + z_cat[0]) * 0.5;
    let expected_cat0 = z_cat[0] - (z_ea[0] + z_b[0]) * 0.5;

    assert!((metrics.ea_imbalance[0] - expected_ea0).abs() < 1e-5);
    assert!((metrics.b_imbalance[0] - expected_b0).abs() < 1e-5);
    assert!((metrics.cat_imbalance[0] - expected_cat0).abs() < 1e-5);
}

#[test]
fn missing_phase_nan() {
    let genesets = vec!["SPLICE_EA_PHASE", "SPLICE_B_PHASE"];
    let values = vec![vec![1.0, 2.0], vec![2.0, 1.0]];
    let activity = make_activity(genesets, values);
    let metrics = compute(&activity).unwrap();

    assert!(metrics.ea_imbalance[0].is_nan());
    assert!(metrics.b_imbalance[0].is_nan());
    assert!(metrics.cat_imbalance[0].is_nan());
}

#[test]
fn deterministic_outputs() {
    let genesets = vec![
        "SPLICE_EA_PHASE",
        "SPLICE_B_PHASE",
        "SPLICE_CATALYTIC_PHASE",
    ];
    let values = vec![vec![1.0, 2.0], vec![2.0, 3.0], vec![3.0, 4.0]];
    let activity = make_activity(genesets, values);

    let first = compute(&activity).unwrap();
    let second = compute(&activity).unwrap();

    assert_eq!(first.ea_imbalance, second.ea_imbalance);
    assert_eq!(first.b_imbalance, second.b_imbalance);
    assert_eq!(first.cat_imbalance, second.cat_imbalance);
    assert_eq!(first.z_ea, second.z_ea);
    assert_eq!(first.z_b, second.z_b);
    assert_eq!(first.z_cat, second.z_cat);
}
