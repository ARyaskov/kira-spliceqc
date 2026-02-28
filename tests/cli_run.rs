use std::fs;
use std::path::Path;

use kira_spliceqc::cli::config::{AnalysisMode, RunConfig, RunMode};
use kira_spliceqc::cli::run::run_pipeline;
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

#[test]
fn end_to_end_json_produced() {
    let input_dir = tempdir().unwrap();
    write_tenx(input_dir.path());

    let out_dir = tempdir().unwrap();
    let config = RunConfig {
        input: input_dir.path().to_path_buf(),
        out_dir: out_dir.path().to_path_buf(),
        cache_path: None,
        mode: AnalysisMode::Cell,
        run_mode: RunMode::Standalone,
        output_json: true,
        output_tsv: false,
        extended: false,
        threads: None,
    };

    run_pipeline(config).unwrap();
    assert!(out_dir.path().join("spliceqc.json").exists());
}

#[test]
fn json_tsv_flag_behavior() {
    let input_dir = tempdir().unwrap();
    write_tenx(input_dir.path());

    let out_dir = tempdir().unwrap();
    let config = RunConfig {
        input: input_dir.path().to_path_buf(),
        out_dir: out_dir.path().to_path_buf(),
        cache_path: None,
        mode: AnalysisMode::Cell,
        run_mode: RunMode::Standalone,
        output_json: true,
        output_tsv: false,
        extended: false,
        threads: None,
    };

    run_pipeline(config).unwrap();
    assert!(out_dir.path().join("spliceqc.json").exists());
    assert!(!out_dir.path().join("spliceqc.tsv").exists());
}

#[test]
fn deterministic_output_hash() {
    let input_dir = tempdir().unwrap();
    write_tenx(input_dir.path());

    let out_dir1 = tempdir().unwrap();
    let out_dir2 = tempdir().unwrap();

    let config1 = RunConfig {
        input: input_dir.path().to_path_buf(),
        out_dir: out_dir1.path().to_path_buf(),
        cache_path: None,
        mode: AnalysisMode::Cell,
        run_mode: RunMode::Standalone,
        output_json: true,
        output_tsv: false,
        extended: false,
        threads: None,
    };
    let config2 = RunConfig {
        input: input_dir.path().to_path_buf(),
        out_dir: out_dir2.path().to_path_buf(),
        cache_path: None,
        mode: AnalysisMode::Cell,
        run_mode: RunMode::Standalone,
        output_json: true,
        output_tsv: false,
        extended: false,
        threads: None,
    };

    run_pipeline(config1).unwrap();
    run_pipeline(config2).unwrap();

    let bytes1 = fs::read(out_dir1.path().join("spliceqc.json")).unwrap();
    let bytes2 = fs::read(out_dir2.path().join("spliceqc.json")).unwrap();

    let mut hasher1 = Sha256::new();
    hasher1.update(&bytes1);
    let hash1 = hasher1.finalize();

    let mut hasher2 = Sha256::new();
    hasher2.update(&bytes2);
    let hash2 = hasher2.finalize();

    assert_eq!(hash1[..], hash2[..]);
}
