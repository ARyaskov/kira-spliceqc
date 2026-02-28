use std::fs;
use std::path::Path;

use kira_spliceqc::cli::config::{AnalysisMode, RunConfig, RunMode};
use kira_spliceqc::cli::run::run_pipeline;
use kira_spliceqc::output::pipeline_contract::spliceqc_header;
use sha2::{Digest, Sha256};
use tempfile::tempdir;

fn write_tenx(dir: &Path) {
    let matrix = "%%MatrixMarket matrix coordinate integer general\n%\n5 2 5\n1 1 1\n2 1 1\n3 2 1\n4 2 1\n5 1 1\n";
    fs::write(dir.join("matrix.mtx"), matrix).unwrap();

    let features = "g1\tSNRPC\ng2\tSF3A1\ng3\tSF3B1\ng4\tSRSF1\ng5\tHNRNPA1\n";
    fs::write(dir.join("features.tsv"), features).unwrap();

    let barcodes = "cell1\ncell2\n";
    fs::write(dir.join("barcodes.tsv"), barcodes).unwrap();
}

fn run_pipeline_contract(input: &Path, out: &Path) {
    let config = RunConfig {
        input: input.to_path_buf(),
        out_dir: out.to_path_buf(),
        cache_path: None,
        mode: AnalysisMode::Cell,
        run_mode: RunMode::Pipeline,
        output_json: false,
        output_tsv: false,
        extended: false,
        threads: None,
    };
    run_pipeline(config).unwrap();
}

#[test]
fn tsv_header_and_column_order() {
    let input = tempdir().unwrap();
    write_tenx(input.path());
    let out = tempdir().unwrap();
    run_pipeline_contract(input.path(), out.path());

    let path = out.path().join("kira-spliceqc").join("spliceqc.tsv");
    let data = fs::read_to_string(path).unwrap();
    let header = data.lines().next().unwrap();
    assert_eq!(header, spliceqc_header());
}

#[test]
fn summary_json_schema() {
    let input = tempdir().unwrap();
    write_tenx(input.path());
    let out = tempdir().unwrap();
    run_pipeline_contract(input.path(), out.path());

    let path = out.path().join("kira-spliceqc").join("summary.json");
    let v: serde_json::Value = serde_json::from_slice(&fs::read(path).unwrap()).unwrap();

    assert_eq!(v["tool"]["name"], "kira-spliceqc");
    assert!(v["tool"]["version"].is_string());
    assert!(v["tool"]["simd"].is_string());
    assert!(v["input"]["n_cells"].is_number());
    assert!(v["input"]["species"].is_string());
    assert!(v["distributions"]["splice_fidelity_index"]["median"].is_number());
    assert!(v["distributions"]["splice_fidelity_index"]["p90"].is_number());
    assert!(v["distributions"]["splice_fidelity_index"]["p99"].is_number());
    assert!(v["distributions"]["stress_splicing_index"]["median"].is_number());
    assert!(v["regimes"]["counts"].is_object());
    assert!(v["regimes"]["fractions"].is_object());
    assert!(v["qc"]["low_confidence_fraction"].is_number());
    assert!(v["qc"]["high_splice_noise_fraction"].is_number());
}

#[test]
fn pipeline_step_json_schema() {
    let input = tempdir().unwrap();
    write_tenx(input.path());
    let out = tempdir().unwrap();
    run_pipeline_contract(input.path(), out.path());

    let path = out.path().join("kira-spliceqc").join("pipeline_step.json");
    let v: serde_json::Value = serde_json::from_slice(&fs::read(path).unwrap()).unwrap();

    assert_eq!(v["tool"]["name"], "kira-spliceqc");
    assert_eq!(v["tool"]["stage"], "splicing");
    assert!(v["tool"]["version"].is_string());
    assert_eq!(v["artifacts"]["summary"], "summary.json");
    assert_eq!(v["artifacts"]["primary_metrics"], "spliceqc.tsv");
    assert_eq!(v["artifacts"]["panels"], "panels_report.tsv");
    assert_eq!(
        v["cell_metrics"],
        serde_json::json!({
            "file": "spliceqc.tsv",
            "id_column": "barcode",
            "regime_column": "regime",
            "confidence_column": "confidence",
            "flag_column": "flags"
        })
    );
    assert_eq!(
        v["regimes"],
        serde_json::json!([
            "HighFidelitySplicing",
            "RegulatedAlternativeSplicing",
            "StressInducedSplicing",
            "SpliceNoiseDominant",
            "SplicingCollapse",
            "Unclassified"
        ])
    );
}

#[test]
fn pipeline_contract_outputs_are_deterministic() {
    let input = tempdir().unwrap();
    write_tenx(input.path());
    let out1 = tempdir().unwrap();
    let out2 = tempdir().unwrap();
    run_pipeline_contract(input.path(), out1.path());
    run_pipeline_contract(input.path(), out2.path());

    let base1 = out1.path().join("kira-spliceqc");
    let base2 = out2.path().join("kira-spliceqc");
    for file in [
        "spliceqc.tsv",
        "summary.json",
        "panels_report.tsv",
        "pipeline_step.json",
    ] {
        let bytes1 = fs::read(base1.join(file)).unwrap();
        let bytes2 = fs::read(base2.join(file)).unwrap();
        let mut h1 = Sha256::new();
        h1.update(&bytes1);
        let mut h2 = Sha256::new();
        h2.update(&bytes2);
        assert_eq!(
            h1.finalize()[..],
            h2.finalize()[..],
            "non-deterministic {file}"
        );
    }
}
