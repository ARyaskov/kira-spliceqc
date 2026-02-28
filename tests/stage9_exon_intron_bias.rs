use kira_spliceqc::model::geneset_activity::GenesetActivityMatrix;
use kira_spliceqc::pipeline::stage9_exon_intron_bias::compute;

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
fn synthetic_bias_sign() {
    let genesets = vec!["SRSF_SR", "HNRNP", "U2AF_AXIS"];
    let values = vec![vec![2.0, 0.0], vec![0.0, 2.0], vec![2.0, 0.0]];
    let activity = make_activity(genesets, values);
    let metrics = compute(&activity).unwrap();

    let z_srsf = robust_z(&activity.values[0..2]);
    let z_hnrnp = robust_z(&activity.values[2..4]);
    let z_u2af = robust_z(&activity.values[4..6]);

    let expected0 = (z_srsf[0] + z_u2af[0]) - z_hnrnp[0];
    let expected1 = (z_srsf[1] + z_u2af[1]) - z_hnrnp[1];

    assert!((metrics.exon_definition_bias[0] - expected0).abs() < 1e-5);
    assert!((metrics.exon_definition_bias[1] - expected1).abs() < 1e-5);
}

#[test]
fn missing_geneset_nan() {
    let genesets = vec!["SRSF_SR", "HNRNP"];
    let values = vec![vec![1.0, 2.0], vec![2.0, 1.0]];
    let activity = make_activity(genesets, values);
    let metrics = compute(&activity).unwrap();

    assert!(metrics.exon_definition_bias[0].is_nan());
    assert!(metrics.exon_definition_bias[1].is_nan());
}

#[test]
fn deterministic_outputs() {
    let genesets = vec!["SRSF_SR", "HNRNP", "U2AF_AXIS"];
    let values = vec![vec![1.0, 2.0], vec![2.0, 3.0], vec![3.0, 4.0]];
    let activity = make_activity(genesets, values);

    let first = compute(&activity).unwrap();
    let second = compute(&activity).unwrap();

    assert_eq!(first.exon_definition_bias, second.exon_definition_bias);
    assert_eq!(first.z_srsf, second.z_srsf);
    assert_eq!(first.z_u2af, second.z_u2af);
    assert_eq!(first.z_hnrnp, second.z_hnrnp);
}
