//! Cell metadata readers (10x metadata.tsv, AnnData obs) and the stratified
//! reference applied to Tier A: deviations are computed within strata and a
//! damaged-cell candidate is flagged only relative to its own cell type.

use std::fs;
use std::path::Path;

use hdf5::types::VarLenUnicode;
use kira_spliceqc::cli::config::{AnalysisMode, RunConfig, RunMode};
use kira_spliceqc::cli::run::run_pipeline;
use kira_spliceqc::input::metadata::{read_metadata_h5ad, read_metadata_tsv};
use tempfile::tempdir;

// ---------------------------------------------------------------------------
// Readers
// ---------------------------------------------------------------------------

#[test]
fn metadata_tsv_is_aligned_to_cell_order() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("metadata.tsv");
    fs::write(
        &path,
        "barcode\tcell_type\tsample\ncellB\tT cell\ts1\ncellA\tB cell\ts1\ncellZ\tNK\ts2\n",
    )
    .unwrap();
    let cells = vec!["cellA".to_string(), "cellB".to_string(), "cellC".to_string()];
    let md = read_metadata_tsv(&path, &cells).unwrap();
    assert_eq!(md.column("cell_type").unwrap(), &["B cell", "T cell", ""]);
    assert_eq!(md.column("sample").unwrap(), &["s1", "s1", ""]);
    assert_eq!(md.cells_without_metadata, 1);
    let (name, _) = md.resolve(&["celltype", "cell_type"]).unwrap();
    assert_eq!(name, "cell_type");

    // Ragged row -> error.
    fs::write(&path, "barcode\tcell_type\ncellA\tT\tx\n").unwrap();
    assert!(read_metadata_tsv(&path, &cells).is_err());
}

fn write_strings(group: &hdf5::Group, name: &str, values: &[&str]) {
    let v: Vec<VarLenUnicode> = values
        .iter()
        .map(|s| unsafe { VarLenUnicode::from_str_unchecked(*s) })
        .collect();
    group
        .new_dataset::<VarLenUnicode>()
        .shape(values.len())
        .create(name)
        .unwrap()
        .write(&v)
        .unwrap();
}

#[test]
fn h5ad_obs_string_and_categorical_columns() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("obs.h5ad");
    let file = hdf5::File::create(&path).unwrap();
    let obs = file.create_group("obs").unwrap();
    write_strings(&obs, "_index", &["c2", "c1", "c3"]);
    write_strings(&obs, "sample", &["s1", "s1", "s2"]);
    // anndata >= 0.8 categorical: codes + categories; -1 = missing.
    let cat = obs.create_group("cell_type").unwrap();
    write_strings(&cat, "categories", &["B", "T"]);
    cat.new_dataset::<i8>()
        .shape(3)
        .create("codes")
        .unwrap()
        .write(&[1i8, 0, -1])
        .unwrap();
    // numeric column: ignored
    obs.new_dataset::<f32>()
        .shape(3)
        .create("n_counts")
        .unwrap()
        .write(&[1.0f32, 2.0, 3.0])
        .unwrap();
    drop(file);

    let obs_order = vec!["c2".to_string(), "c1".to_string(), "c3".to_string()];
    let cells = vec!["c1".to_string(), "c2".to_string(), "c3".to_string()];
    let md = read_metadata_h5ad(&path, &obs_order, &cells).unwrap();
    assert_eq!(md.column("sample").unwrap(), &["s1", "s1", "s2"]);
    assert_eq!(md.column("cell_type").unwrap(), &["B", "T", ""]);
    assert!(md.column("n_counts").is_none());
    assert!(md.source.starts_with("h5ad-obs:"));
}

// ---------------------------------------------------------------------------
// Stratified Tier A end-to-end
// ---------------------------------------------------------------------------

/// 5 panel genes x 120 cells: 60 "Neuron" cells with UF ~ 0.6, 60 "Glia"
/// cells with UF ~ 0.2. One glia cell is damaged (UF 0.02) and one neuron
/// cell has UF 0.2, which is normal for glia but 4 sigma low for neurons.
fn write_stratified_dataset(dir: &Path) {
    fs::create_dir_all(dir).unwrap();
    let genes = ["SNRPC", "SF3A1", "SF3B1", "SRSF1", "HNRNPA1"];
    let n_cells = 120;
    let mut mtx = String::new();
    let mut spliced = String::new();
    let mut unspliced = String::new();
    let mut barcodes = String::new();
    let mut metadata = String::from("barcode\tcell_type\n");
    let mut entries = 0;
    let mut body_m = String::new();
    let mut body_s = String::new();
    let mut body_u = String::new();
    for c in 0..n_cells {
        let neuron = c < 60;
        let barcode = format!("CELL{c:03}");
        barcodes.push_str(&format!("{barcode}\n"));
        metadata.push_str(&format!("{barcode}\t{}\n", if neuron { "Neuron" } else { "Glia" }));
        // per-gene counts: 200 UMIs per gene with a per-cell jitter so MAD > 0
        let jitter = (c % 5) as f64 * 0.01;
        let uf = match c {
            0 => 0.2,          // neuron with glia-like UF: outlier in its stratum
            60 => 0.02,        // damaged glia cell
            _ if neuron => 0.6 + jitter,
            _ => 0.2 + jitter,
        };
        for g in 1..=genes.len() {
            let u = (200.0 * uf).round() as u32;
            // Per-gene jitter on the spliced count so that panel scores are
            // not identical across cells (a zero MAD gives undefined z-scores).
            let s = 200 - u + ((c * 7 + g) % 5) as u32;
            body_m.push_str(&format!("{g} {} {}\n", c + 1, s + u));
            body_s.push_str(&format!("{g} {} {s}\n", c + 1));
            body_u.push_str(&format!("{g} {} {u}\n", c + 1));
            entries += 1;
        }
    }
    let header = format!(
        "%%MatrixMarket matrix coordinate integer general\n{} {} {}\n",
        genes.len(),
        n_cells,
        entries
    );
    mtx.push_str(&header);
    mtx.push_str(&body_m);
    spliced.push_str(&header);
    spliced.push_str(&body_s);
    unspliced.push_str(&header);
    unspliced.push_str(&body_u);
    fs::write(dir.join("matrix.mtx"), mtx).unwrap();
    fs::write(dir.join("spliced.mtx"), spliced).unwrap();
    fs::write(dir.join("unspliced.mtx"), unspliced).unwrap();
    fs::write(
        dir.join("features.tsv"),
        genes.iter().enumerate().map(|(i, g)| format!("g{i}\t{g}\n")).collect::<String>(),
    )
    .unwrap();
    fs::write(dir.join("barcodes.tsv"), barcodes).unwrap();
    fs::write(dir.join("metadata.tsv"), metadata).unwrap();
}

fn run(input: &Path, out: &Path, stratify_by: Option<&str>) {
    run_pipeline(RunConfig {
        input: input.to_path_buf(),
        out_dir: out.to_path_buf(),
        cache_path: None,
        layers: None,
        metadata: None,
        stratify_by: stratify_by.map(str::to_string),
        reference: None,
        mode: AnalysisMode::Cell,
        run_mode: RunMode::Pipeline,
        output_json: true,
        output_tsv: true,
        extended: false,
        threads: None,
        experimental_signatures: false,
    })
    .unwrap();
}

fn read_cells(path: &Path) -> (Vec<String>, Vec<Vec<String>>) {
    let text = fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let header = lines.next().unwrap().split('\t').map(str::to_string).collect();
    let rows = lines
        .map(|l| l.split('\t').map(str::to_string).collect())
        .collect();
    (header, rows)
}

#[test]
fn stratified_reference_flags_outliers_within_cell_type() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("data");
    write_stratified_dataset(&input);
    let out = tempdir().unwrap();
    run(&input, out.path(), None);
    let base = out.path().join("kira-spliceqc");

    let summary: serde_json::Value =
        serde_json::from_slice(&fs::read(base.join("summary.json")).unwrap()).unwrap();
    assert_eq!(summary["reference"]["mode"], "stratified");
    assert_eq!(summary["reference"]["column"], "cell_type");
    assert_eq!(summary["reference"]["n_strata"], 3); // global + Neuron + Glia
    assert_eq!(summary["reference"]["folded_cells"], 0);
    let strata = summary["unspliced"]["strata"].as_array().unwrap();
    let glia = strata.iter().find(|s| s["name"] == "Glia").unwrap();
    let neuron = strata.iter().find(|s| s["name"] == "Neuron").unwrap();
    assert!((glia["median"].as_f64().unwrap() - 0.22).abs() < 0.03);
    assert!((neuron["median"].as_f64().unwrap() - 0.62).abs() < 0.03);

    let (header, rows) = read_cells(&base.join("cells.tsv"));
    let col = |n: &str| header.iter().position(|h| h == n).unwrap();
    let (name, dev, flag) = (
        col("cell_name"),
        col("unspliced_fraction_dev"),
        col("nuclear_fraction_flag"),
    );
    let mut flagged = Vec::new();
    for r in &rows {
        if r[flag] == "true" {
            flagged.push(r[name].clone());
        }
        match r[name].as_str() {
            "CELL000" => assert!(r[dev].parse::<f64>().unwrap() < -3.0, "neuron outlier {}", r[dev]),
            "CELL060" => assert!(r[dev].parse::<f64>().unwrap() < -3.0, "damaged glia {}", r[dev]),
            "CELL061" => assert!(r[dev].parse::<f64>().unwrap().abs() < 2.0, "typical glia {}", r[dev]),
            _ => {}
        }
    }
    flagged.sort();
    assert_eq!(flagged, vec!["CELL000", "CELL060"]);

    let cells: serde_json::Value =
        serde_json::from_slice(&fs::read(base.join("cells.json")).unwrap()).unwrap();
    assert_eq!(cells["reference"]["mode"], "stratified");
    assert_eq!(cells["reference"]["labels"].as_array().unwrap().len(), 120);
}

#[test]
fn global_reference_on_mixed_cell_types_misses_outliers() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("data");
    write_stratified_dataset(&input);
    fs::remove_file(input.join("metadata.tsv")).unwrap();
    let out = tempdir().unwrap();
    run(&input, out.path(), None);
    let base = out.path().join("kira-spliceqc");

    let summary: serde_json::Value =
        serde_json::from_slice(&fs::read(base.join("summary.json")).unwrap()).unwrap();
    assert_eq!(summary["reference"]["mode"], "global");
    assert!(summary["reference"]["column"].is_null());

    // One global stratum over a bimodal population (UF 0.6 vs 0.2): the
    // between-type spread inflates the overdispersion so much that neither
    // the neuron at UF 0.2 nor even the damaged cell at 0.02 is an outlier.
    // This is why the reference must be stratified (METRICS.md).
    let (header, rows) = read_cells(&base.join("cells.tsv"));
    let col = |n: &str| header.iter().position(|h| h == n).unwrap();
    let flagged: Vec<&str> = rows
        .iter()
        .filter(|r| r[col("nuclear_fraction_flag")] == "true")
        .map(|r| r[col("cell_name")].as_str())
        .collect();
    assert!(flagged.is_empty(), "{flagged:?}");
    let damaged = rows.iter().find(|r| r[col("cell_name")] == "CELL060").unwrap();
    let dev: f64 = damaged[col("unspliced_fraction_dev")].parse().unwrap();
    assert!(dev > -3.0 && dev < 0.0, "{dev}");
}

#[test]
fn unknown_stratify_column_falls_back_to_global() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("data");
    write_stratified_dataset(&input);
    let out = tempdir().unwrap();
    run(&input, out.path(), Some("does_not_exist"));
    let summary: serde_json::Value = serde_json::from_slice(
        &fs::read(out.path().join("kira-spliceqc").join("summary.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(summary["reference"]["mode"], "global");
}
