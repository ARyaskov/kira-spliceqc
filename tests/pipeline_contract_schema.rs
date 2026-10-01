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
        layers: None,
        junctions: None,
        metadata: None,
        stratify_by: None,
        reference: None,
        catalog: None,
        min_counts: 0,
        min_genes: 0,
        mode: AnalysisMode::Cell,
        run_mode: RunMode::Pipeline,
        output_json: false,
        output_tsv: false,
        extended: false,
        threads: None,
        experimental_signatures: false,
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
    // Distribution statistics are `null` when no cell has a finite value
    // (missing data is never coerced to 0).
    for metric in ["splice_fidelity_index", "stress_splicing_index"] {
        for stat in ["median", "p90", "p99"] {
            let value = &v["distributions"][metric][stat];
            assert!(value.is_number() || value.is_null(), "{metric}.{stat}");
        }
    }
    assert!(v["regimes"]["counts"].is_object());
    assert!(v["regimes"]["fractions"].is_object());
    assert!(v["qc"]["low_confidence_fraction"].is_number());
    assert!(v["qc"]["high_splice_noise_fraction"].is_number());
    assert_eq!(
        v["splicing_instability"]["panel_version"],
        "SPLICEQC_INSTABILITY_PANEL_V1"
    );
    assert!(v["splicing_instability"]["thresholds"]["flag_rule"].is_string());
    assert_eq!(
        v["splicing_instability"]["thresholds"]["deviation_threshold"],
        3.0
    );
    assert!(
        v["splicing_instability"]["global_stats"]["sos_p50"].is_number()
            || v["splicing_instability"]["global_stats"]["sos_p50"].is_null()
    );
    assert!(v["splicing_instability"]["missingness"]["panel_coverage"].is_array());

    // Provenance: tool, command, catalog hash, parameters, undefined counts.
    let p = &v["provenance"];
    assert_eq!(p["tool"]["name"], "kira-spliceqc");
    assert_eq!(p["command"]["run_mode"], "pipeline");
    assert_eq!(p["command"]["experimental_signatures"], true);
    assert!(p["geneset_catalog"]["crc64"].as_str().unwrap().len() == 16);
    assert!(p["geneset_catalog"]["source"].is_string());
    assert_eq!(p["input_levels"], serde_json::json!(["L0"]));
    assert_eq!(p["reference"]["mode"], "global");
    assert_eq!(p["parameters"]["controls_per_gene"], 50);
    assert_eq!(p["parameters"]["min_stratum_cells"], 50);
    assert!(p["undefined_cells"]["sis"].is_number());
    assert!(p["undefined_cells"]["unspliced_fraction"].is_null());
    assert!(p.get("reference_file").is_none());
}

#[test]
fn multiqc_custom_content() {
    let input = tempdir().unwrap();
    write_tenx(input.path());
    let out = tempdir().unwrap();
    run_pipeline_contract(input.path(), out.path());

    let path = out
        .path()
        .join("kira-spliceqc")
        .join("kira_spliceqc_mqc.json");
    let v: serde_json::Value = serde_json::from_slice(&fs::read(path).unwrap()).unwrap();
    assert_eq!(v["id"], "kira_spliceqc");
    assert_eq!(v["plot_type"], "table");
    let sample = input.path().file_name().unwrap().to_str().unwrap();
    let row = &v["data"][sample];
    assert_eq!(row["n_cells"], 2);
    assert!(row["low_depth_fraction"].is_number());
    assert_eq!(row["reference_mode"], "global");
    // No layers, no junctions: no Tier A / Tier B columns.
    assert!(row.get("unspliced_fraction_median").is_none());
    assert!(row.get("cryptic_3ss_high_fraction").is_none());
    assert!(v["headers"]["n_cells"]["title"].is_string());
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
    assert_eq!(v["artifacts"]["multiqc"], "kira_spliceqc_mqc.json");
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
        "kira_spliceqc_mqc.json",
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
