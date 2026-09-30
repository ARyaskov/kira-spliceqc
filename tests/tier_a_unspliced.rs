//! End-to-end Tier A: a 10x directory with spliced/unspliced layers produces
//! per-cell unspliced fractions in cells.tsv / cells.json and a Tier A block
//! in the pipeline summary.

use std::fs;
use std::path::Path;

use kira_spliceqc::cli::config::{AnalysisMode, RunConfig, RunMode};
use kira_spliceqc::cli::run::run_pipeline;
use tempfile::tempdir;

/// 5 genes (all in the splicing panels so the L0 stages run) x 3 cells.
fn write_dataset(dir: &Path, with_layers: bool) {
    fs::create_dir_all(dir).unwrap();
    let genes = ["SNRPC", "SF3A1", "SF3B1", "SRSF1", "HNRNPA1"];
    let mut mtx = String::from("%%MatrixMarket matrix coordinate integer general\n5 3 15\n");
    let mut spliced = String::from("%%MatrixMarket matrix coordinate integer general\n5 3 15\n");
    let mut unspliced = String::from("%%MatrixMarket matrix coordinate integer general\n5 3 15\n");
    // cell 1: ~300 spliced + ~100 unspliced per gene (UF exactly 0.25)
    // cell 2: ~100 spliced + ~100 unspliced per gene (UF exactly 0.5)
    // cell 3: 10 spliced + 5 unspliced per gene  (< 100 UMIs -> undefined)
    // Per-gene jitter keeps cell profiles non-proportional (identical cp10k
    // profiles give zero-MAD panels and undefined z-scores) while preserving
    // each cell's overall unspliced fraction.
    for g in 1..=5u32 {
        for c in 0..3u32 {
            // jitter differs per gene *and* per cell
            let j = (g * 7 + c * 3) % 5;
            let (s, u) = match c {
                0 => (300 + 3 * j, 100 + j),
                1 => (100 + 2 * j, 100 + 2 * j),
                _ => (10 + j, 5),
            };
            mtx.push_str(&format!("{g} {} {}\n", c + 1, s + u));
            spliced.push_str(&format!("{g} {} {s}\n", c + 1));
            unspliced.push_str(&format!("{g} {} {u}\n", c + 1));
        }
    }
    fs::write(dir.join("matrix.mtx"), mtx).unwrap();
    let features = genes
        .iter()
        .enumerate()
        .map(|(i, g)| format!("g{i}\t{g}\n"))
        .collect::<String>();
    fs::write(dir.join("features.tsv"), features).unwrap();
    fs::write(dir.join("barcodes.tsv"), "cellA\ncellB\ncellC\n").unwrap();
    if with_layers {
        fs::write(dir.join("spliced.mtx"), spliced).unwrap();
        fs::write(dir.join("unspliced.mtx"), unspliced).unwrap();
    }
}

fn config(input: &Path, out: &Path, run_mode: RunMode) -> RunConfig {
    RunConfig {
        input: input.to_path_buf(),
        out_dir: out.to_path_buf(),
        cache_path: None,
        layers: None,
        metadata: None,
        stratify_by: None,
        reference: None,
        mode: AnalysisMode::Cell,
        run_mode,
        output_json: true,
        output_tsv: true,
        extended: false,
        threads: None,
        experimental_signatures: false,
    }
}

fn column(header: &str, name: &str) -> usize {
    header.split('\t').position(|h| h == name).unwrap_or_else(|| panic!("missing column {name}"))
}

#[test]
fn unspliced_fraction_is_written_per_cell() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("data");
    write_dataset(&input, true);
    let out = tempdir().unwrap();
    run_pipeline(config(&input, out.path(), RunMode::Standalone)).unwrap();

    let tsv = fs::read_to_string(out.path().join("cells.tsv")).unwrap();
    let mut lines = tsv.lines();
    let header = lines.next().unwrap();
    let uf = column(header, "unspliced_fraction");
    let lo = column(header, "unspliced_fraction_ci_low");
    let hi = column(header, "unspliced_fraction_ci_high");
    let su = column(header, "spliced_umis");
    let name = column(header, "cell_name");
    // Composite columns are absent without --experimental-signatures.
    assert!(!header.split('\t').any(|h| h == "sis" || h == "SOS"));

    let mut seen = 0;
    for line in lines {
        let f: Vec<&str> = line.split('\t').collect();
        match f[name] {
            "cellA" => {
                // 5 genes x (300 + 3j), j = 7g mod 5 = 2, 4, 1, 3, 0 -> 1530
                assert_eq!(f[su], "1530");
                assert!((f[uf].parse::<f64>().unwrap() - 0.25).abs() < 1e-6);
                assert!(f[lo].parse::<f64>().unwrap() < 0.25);
                assert!(f[hi].parse::<f64>().unwrap() > 0.25);
            }
            "cellB" => assert!((f[uf].parse::<f64>().unwrap() - 0.5).abs() < 1e-6),
            "cellC" => {
                assert_eq!(f[su], "60");
                assert_eq!(f[uf], "", "75 layer UMIs must be undefined");
            }
            other => panic!("unexpected cell {other}"),
        }
        seen += 1;
    }
    assert_eq!(seen, 3);

    let v: serde_json::Value =
        serde_json::from_slice(&fs::read(out.path().join("cells.json")).unwrap()).unwrap();
    assert_eq!(v["input_levels"], serde_json::json!(["L0", "L1"]));
    assert_eq!(v["unspliced"]["undefined_cells"], 1);
    // 3 cells < MIN_CELLS_PER_GENE (10): no gene has a reference ratio, so
    // the intron retention index is undefined; the block is still present.
    assert_eq!(v["intron_retention"]["undefined_cells"], 3);
    assert_eq!(v["intron_retention"]["genes_with_reference"], 0);
    assert_eq!(v["cells"][0]["intron_retention"]["ir_genes_used"], 0);
    let iri = column(header, "intron_retention_index");
    assert!(tsv.lines().skip(1).all(|l| l.split('\t').nth(iri).unwrap().is_empty()));
    assert!(v["unspliced"]["source"].as_str().unwrap().starts_with("mtx-dir:"));
    assert!(v["cells"][0]["unspliced"]["spliced_umis"].is_number());
}

#[test]
fn without_layers_tier_a_columns_are_empty_and_summary_says_l0() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("data");
    write_dataset(&input, false);
    let out = tempdir().unwrap();
    run_pipeline(config(&input, out.path(), RunMode::Pipeline)).unwrap();

    let base = out.path().join("kira-spliceqc");
    let tsv = fs::read_to_string(base.join("cells.tsv")).unwrap();
    let mut lines = tsv.lines();
    let header = lines.next().unwrap();
    let uf = column(header, "unspliced_fraction");
    for line in lines {
        assert_eq!(line.split('\t').nth(uf).unwrap(), "");
    }
    let v: serde_json::Value =
        serde_json::from_slice(&fs::read(base.join("summary.json")).unwrap()).unwrap();
    assert_eq!(v["input"]["levels"], serde_json::json!(["L0"]));
    assert!(v["unspliced"].is_null());
    let cells: serde_json::Value =
        serde_json::from_slice(&fs::read(base.join("cells.json")).unwrap()).unwrap();
    assert!(cells["unspliced"].is_null());
    assert!(cells["cells"][0].get("unspliced").is_none());
}

#[test]
fn pipeline_summary_carries_tier_a_block() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("data");
    write_dataset(&input, true);
    let out = tempdir().unwrap();
    run_pipeline(config(&input, out.path(), RunMode::Pipeline)).unwrap();

    let v: serde_json::Value = serde_json::from_slice(
        &fs::read(out.path().join("kira-spliceqc").join("summary.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(v["input"]["levels"], serde_json::json!(["L0", "L1"]));
    assert_eq!(v["unspliced"]["n_defined_cells"], 2);
    assert_eq!(v["unspliced"]["min_layer_umis"], 100);
    // Two defined cells (0.25, 0.5); the round((n-1)q) quantile rule picks
    // one of them, so the median is one of the two values.
    let median = v["unspliced"]["median"].as_f64().unwrap();
    assert!((median - 0.25).abs() < 1e-6 || (median - 0.5).abs() < 1e-6, "{median}");
    assert!((v["unspliced"]["p10"].as_f64().unwrap() - 0.25).abs() < 1e-6);
    assert!((v["unspliced"]["p90"].as_f64().unwrap() - 0.5).abs() < 1e-6);
}
