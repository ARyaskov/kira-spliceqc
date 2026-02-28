use kira_spliceqc::model::geneset_activity::GenesetActivityMatrix;
use kira_spliceqc::pipeline::stage11_splicing_noise::compute;

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

fn noise_for(
    metrics: &kira_spliceqc::model::splicing_noise::SplicingNoiseMetrics,
    id: &str,
) -> f32 {
    metrics
        .per_geneset_noise
        .iter()
        .find(|(name, _)| name == id)
        .map(|(_, v)| *v)
        .unwrap_or(f32::NAN)
}

#[test]
fn stable_dataset_low_noise() {
    let genesets = vec!["U1_CORE", "U2_CORE"];
    let values = vec![vec![2.0, 2.0, 2.0, 2.0], vec![1.0, 1.0, 1.0, 1.0]];
    let activity = make_activity(genesets, values);
    let metrics = compute(&activity).unwrap();

    assert!(metrics.noise_index.abs() < 1e-6);
    assert!(noise_for(&metrics, "U1_CORE").abs() < 1e-6);
    assert!(noise_for(&metrics, "U2_CORE").abs() < 1e-6);
}

#[test]
fn noisy_dataset_higher_noise() {
    let genesets = vec!["U1_CORE", "U2_CORE"];
    let values = vec![vec![2.0, 2.0, 2.0, 2.0], vec![0.0, 10.0, 0.0, 10.0]];
    let activity = make_activity(genesets, values);
    let metrics = compute(&activity).unwrap();

    let stable = noise_for(&metrics, "U1_CORE");
    let noisy = noise_for(&metrics, "U2_CORE");

    assert!(noisy > stable);
    assert!(noisy > 0.5);
    assert!(metrics.noise_index > 0.2);
}

#[test]
fn missing_genesets_exclude_nan() {
    let genesets = vec!["U1_CORE", "U2_CORE"];
    let values = vec![vec![0.0, 0.0, 0.0], vec![1.0, 1.0, 1.0]];
    let activity = make_activity(genesets, values);
    let metrics = compute(&activity).unwrap();

    assert!(noise_for(&metrics, "U1_CORE").is_nan());
    assert!(metrics.noise_index.abs() < 1e-6);
}

#[test]
fn deterministic_outputs() {
    let genesets = vec!["U1_CORE", "U2_CORE", "SF3B_AXIS"];
    let values = vec![
        vec![1.0, 2.0, 3.0, 4.0],
        vec![2.0, 3.0, 4.0, 5.0],
        vec![3.0, 4.0, 5.0, 6.0],
    ];
    let activity = make_activity(genesets, values);

    let first = compute(&activity).unwrap();
    let second = compute(&activity).unwrap();

    assert_eq!(first.noise_index, second.noise_index);
    assert_eq!(first.per_geneset_noise, second.per_geneset_noise);
}
