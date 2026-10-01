//! Cell QC flags: LOW_DEPTH from --min-counts / --min-genes and DOUBLET from
//! a metadata column; both are excluded from reference norms and get
//! undefined deviations.

use std::fs;
use std::path::Path;

use kira_spliceqc::cli::config::{AnalysisMode, RunConfig, RunMode};
use kira_spliceqc::cli::run::run_pipeline;
use tempfile::tempdir;

/// 25 genes x 100 cells with layers; cells 0..10 are shallow (~60 UMIs),
/// cells 10..15 are marked doublets in metadata.tsv.
fn write_dataset(dir: &Path) {
    fs::create_dir_all(dir).unwrap();
    let mut genes: Vec<String> = ["SNRPC", "SF3A1", "SF3B1", "SRSF1", "HNRNPA1"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    genes.extend((0..20).map(|i| format!("FILLER{i:02}")));
    let n_cells = 100;
    let (mut m, mut s, mut u, mut barcodes) = (String::new(), String::new(), String::new(), String::new());
    let mut metadata = String::from("barcode\tpredicted_doublet\n");
    let mut entries = 0;
    for c in 0..n_cells {
        let barcode = format!("CELL{c:03}");
        barcodes.push_str(&format!("{barcode}\n"));
        metadata.push_str(&format!("{barcode}\t{}\n", if (10..15).contains(&c) { "True" } else { "False" }));
        for g in 1..=genes.len() {
            let total = if c < 10 { 2 + (g % 3) as u32 } else { 200 + ((c * 11 + g * 7) % 9) as u32 };
            let un = (total as f64 * (0.25 + (c % 5) as f64 * 0.01)).round() as u32;
            m.push_str(&format!("{g} {} {total}\n", c + 1));
            s.push_str(&format!("{g} {} {}\n", c + 1, total - un));
            u.push_str(&format!("{g} {} {un}\n", c + 1));
            entries += 1;
        }
    }
    let header = format!("%%MatrixMarket matrix coordinate integer general\n{} {} {}\n", genes.len(), n_cells, entries);
    fs::write(dir.join("matrix.mtx"), format!("{header}{m}")).unwrap();
    fs::write(dir.join("spliced.mtx"), format!("{header}{s}")).unwrap();
    fs::write(dir.join("unspliced.mtx"), format!("{header}{u}")).unwrap();
    fs::write(
        dir.join("features.tsv"),
        genes.iter().enumerate().map(|(i, g)| format!("g{i}\t{g}\n")).collect::<String>(),
    )
    .unwrap();
    fs::write(dir.join("barcodes.tsv"), barcodes).unwrap();
    fs::write(dir.join("metadata.tsv"), metadata).unwrap();
}

#[test]
fn low_depth_and_doublet_cells_are_flagged_and_excluded() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("data");
    write_dataset(&input);
    let out = tempdir().unwrap();
    run_pipeline(RunConfig {
        input: input.clone(),
        out_dir: out.path().to_path_buf(),
        cache_path: None,
        layers: None,
        junctions: None,
        metadata: None,
        stratify_by: None,
        reference: None,
        catalog: None,
        min_counts: 500,
        min_genes: 20,
        mode: AnalysisMode::Cell,
        run_mode: RunMode::Pipeline,
        output_json: true,
        output_tsv: true,
        extended: false,
        threads: None,
        experimental_signatures: false,
    })
    .unwrap();
    let base = out.path().join("kira-spliceqc");

    let cells = fs::read_to_string(base.join("cells.tsv")).unwrap();
    let mut lines = cells.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let col = |n: &str| header.iter().position(|h| *h == n).unwrap();
    let (name, low, dbl, uf_dev, iri_dev) = (
        col("cell_name"),
        col("low_depth"),
        col("doublet"),
        col("unspliced_fraction_dev"),
        col("intron_retention_index_dev"),
    );
    let mut n_low = 0;
    let mut n_dbl = 0;
    for line in lines {
        let f: Vec<&str> = line.split('\t').collect();
        let idx: usize = f[name][4..].parse().unwrap();
        if f[low] == "true" {
            n_low += 1;
            assert!(idx < 10, "{}", f[name]);
            assert_eq!(f[uf_dev], "", "excluded cells have no deviation");
            assert_eq!(f[iri_dev], "");
        }
        if f[dbl] == "true" {
            n_dbl += 1;
            assert!((10..15).contains(&idx), "{}", f[name]);
            assert_eq!(f[uf_dev], "");
        }
        if idx >= 15 {
            assert!(!f[uf_dev].is_empty(), "kept cell {} must have a deviation", f[name]);
        }
    }
    assert_eq!(n_low, 10);
    assert_eq!(n_dbl, 5);

    let contract = fs::read_to_string(base.join("spliceqc.tsv")).unwrap();
    assert!(contract.lines().any(|l| l.starts_with("CELL000\t") && l.contains("LOW_DEPTH")));
    assert!(contract.lines().any(|l| l.starts_with("CELL012\t") && l.contains("DOUBLET")));

    let summary: serde_json::Value =
        serde_json::from_slice(&fs::read(base.join("summary.json")).unwrap()).unwrap();
    assert!((summary["qc"]["low_depth_fraction"].as_f64().unwrap() - 0.10).abs() < 1e-9);
    assert!((summary["qc"]["doublet_fraction"].as_f64().unwrap() - 0.05).abs() < 1e-9);
    assert_eq!(summary["qc"]["doublet_column"], "predicted_doublet");
    assert_eq!(summary["provenance"]["reference"]["excluded_cells"], 15);
    assert_eq!(summary["provenance"]["parameters"]["min_counts"], 500);
    // Excluded cells are not part of any stratum's norm.
    assert_eq!(summary["reference"]["strata"][0]["n_cells"], 85);
}
