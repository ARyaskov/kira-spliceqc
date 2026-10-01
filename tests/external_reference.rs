//! External reference: `reference build` on a control dataset, then a run of
//! a second dataset with `--reference`. A whole stratum that shifted relative
//! to the control is invisible to an internal (dataset-relative) reference
//! and obvious under the external one.

use std::fs;
use std::path::Path;

use kira_spliceqc::cli::config::{AnalysisMode, RunConfig, RunMode};
use kira_spliceqc::cli::run::{build_reference_file, run_pipeline};
use kira_spliceqc::reference::external::ReferenceFile;
use tempfile::tempdir;

/// 25 genes (5 splicing-panel genes + 20 fillers) x 120 cells: 60 "Neuron"
/// cells at unspliced fraction 0.6 and 60 "Glia" cells at `glia_uf`.
fn write_dataset(dir: &Path, glia_uf: f64) {
    fs::create_dir_all(dir).unwrap();
    let mut genes: Vec<String> = ["SNRPC", "SF3A1", "SF3B1", "SRSF1", "HNRNPA1"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    genes.extend((0..20).map(|i| format!("FILLER{i:02}")));
    let n_cells = 120;
    let (mut m, mut s, mut u, mut barcodes, mut metadata) =
        (String::new(), String::new(), String::new(), String::new(), String::from("barcode\tcell_type\n"));
    let mut entries = 0;
    for c in 0..n_cells {
        let neuron = c < 60;
        let barcode = format!("CELL{c:03}");
        barcodes.push_str(&format!("{barcode}\n"));
        metadata.push_str(&format!("{barcode}\t{}\n", if neuron { "Neuron" } else { "Glia" }));
        let uf = if neuron { 0.6 } else { glia_uf } + (c % 5) as f64 * 0.01;
        for g in 1..=genes.len() {
            let total = 200 + ((c * 11 + g * 7) % 9) as u32;
            let un = (total as f64 * uf).round() as u32;
            let sp = total - un;
            m.push_str(&format!("{g} {} {total}\n", c + 1));
            s.push_str(&format!("{g} {} {sp}\n", c + 1));
            u.push_str(&format!("{g} {} {un}\n", c + 1));
            entries += 1;
        }
    }
    let header = format!(
        "%%MatrixMarket matrix coordinate integer general\n{} {} {}\n",
        genes.len(),
        n_cells,
        entries
    );
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

fn config(input: &Path, out: &Path, reference: Option<&Path>) -> RunConfig {
    RunConfig {
        input: input.to_path_buf(),
        out_dir: out.to_path_buf(),
        cache_path: None,
        layers: None,
        junctions: None,
        metadata: None,
        stratify_by: None,
        reference: reference.map(Path::to_path_buf),
        catalog: None,
        min_counts: 0,
        min_genes: 0,
        mode: AnalysisMode::Cell,
        run_mode: RunMode::Pipeline,
        output_json: true,
        output_tsv: true,
        extended: false,
        threads: None,
        experimental_signatures: false,
    }
}

struct Cells {
    header: Vec<String>,
    rows: Vec<Vec<String>>,
}

impl Cells {
    fn read(path: &Path) -> Self {
        let text = fs::read_to_string(path).unwrap();
        let mut lines = text.lines();
        let header = lines.next().unwrap().split('\t').map(str::to_string).collect();
        let rows = lines.map(|l| l.split('\t').map(str::to_string).collect()).collect();
        Self { header, rows }
    }
    fn median(&self, col: &str, cells: impl Fn(&str) -> bool) -> f64 {
        let ci = self.header.iter().position(|h| h == col).unwrap();
        let ni = self.header.iter().position(|h| h == "cell_name").unwrap();
        let mut v: Vec<f64> = self
            .rows
            .iter()
            .filter(|r| cells(&r[ni]))
            .filter_map(|r| r[ci].parse().ok())
            .collect();
        assert!(!v.is_empty(), "no defined {col}");
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v[v.len() / 2]
    }
    fn count(&self, col: &str, value: &str, cells: impl Fn(&str) -> bool) -> usize {
        let ci = self.header.iter().position(|h| h == col).unwrap();
        let ni = self.header.iter().position(|h| h == "cell_name").unwrap();
        self.rows.iter().filter(|r| cells(&r[ni]) && r[ci] == value).count()
    }
}

fn glia(name: &str) -> bool {
    name[4..].parse::<usize>().unwrap() >= 60
}
fn neuron(name: &str) -> bool {
    !glia(name)
}

#[test]
fn reference_build_then_apply_reveals_a_shifted_stratum() {
    let dir = tempdir().unwrap();
    let control = dir.path().join("control");
    write_dataset(&control, 0.2);
    let ref_path = dir.path().join("ref.json");
    build_reference_file(config(&control, dir.path(), None), &ref_path).unwrap();

    let file = ReferenceFile::read(&ref_path).unwrap();
    assert_eq!(file.stratify_by.as_deref(), Some("cell_type"));
    assert_eq!(file.strata.len(), 3);
    let glia_ref = file.strata.iter().find(|s| s.name == "Glia").unwrap();
    assert!((glia_ref.unspliced_fraction.unwrap().median - 0.22).abs() < 0.03);
    assert_eq!(glia_ref.gene_unspliced_ratio.len(), 25);
    assert!(glia_ref.intron_retention_index.is_some());

    // Target: glia shifted to UF 0.4 (every gene), neurons unchanged.
    let target = dir.path().join("target");
    write_dataset(&target, 0.4);

    let out_ext = tempdir().unwrap();
    run_pipeline(config(&target, out_ext.path(), Some(&ref_path))).unwrap();
    let base = out_ext.path().join("kira-spliceqc");
    let summary: serde_json::Value =
        serde_json::from_slice(&fs::read(base.join("summary.json")).unwrap()).unwrap();
    assert_eq!(summary["reference"]["mode"], "external");
    assert_eq!(summary["reference"]["column"], "cell_type");
    assert!(summary["reference"]["external_file"].as_str().unwrap().ends_with("ref.json"));
    assert_eq!(
        summary["reference"]["external_metrics"],
        serde_json::json!(["unspliced_fraction", "intron_retention_index"])
    );
    let ext = Cells::read(&base.join("cells.tsv"));
    let glia_uf_dev = ext.median("unspliced_fraction_dev", glia);
    let neuron_uf_dev = ext.median("unspliced_fraction_dev", neuron);
    assert!(glia_uf_dev > 5.0, "glia UF deviation under external reference: {glia_uf_dev}");
    assert!(neuron_uf_dev.abs() < 1.5, "neuron UF deviation: {neuron_uf_dev}");
    // Intron retention doubled in glia relative to the control: index ~ log2(2).
    let glia_iri = ext.median("intron_retention_index", glia);
    assert!(glia_iri > 0.6 && glia_iri < 1.2, "glia IRI vs control: {glia_iri}");
    assert!(ext.median("intron_retention_index", neuron).abs() < 0.2);
    assert!(ext.median("intron_retention_index_dev", glia) > 3.0);
    assert!(ext.count("intron_retention_high", "true", glia) >= 50);
    assert_eq!(ext.count("nuclear_fraction_flag", "true", glia), 0);

    // Same target with the internal (dataset-relative) reference: the shift
    // is the new normal and nothing deviates.
    let out_int = tempdir().unwrap();
    run_pipeline(config(&target, out_int.path(), None)).unwrap();
    let int = Cells::read(&out_int.path().join("kira-spliceqc").join("cells.tsv"));
    assert!(int.median("unspliced_fraction_dev", glia).abs() < 1.0);
    assert!(int.median("intron_retention_index", glia).abs() < 0.2);
    assert_eq!(int.count("intron_retention_high", "true", glia), 0);
}

#[test]
fn reference_file_is_validated() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("bad.json");
    fs::write(&path, r#"{"format":"something-else","version":1,"tool_version":"x","stratify_by":null,"n_cells":1,"min_layer_umis":100,"strata":[]}"#).unwrap();
    let err = ReferenceFile::read(&path).unwrap_err().to_string();
    assert!(err.contains("format"), "{err}");
    fs::write(&path, r#"{"format":"kira-spliceqc-reference","version":1,"tool_version":"x","stratify_by":null,"n_cells":1,"min_layer_umis":100,"strata":[]}"#).unwrap();
    let err = ReferenceFile::read(&path).unwrap_err().to_string();
    assert!(err.contains("global"), "{err}");
}

#[test]
fn reference_build_requires_layers() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("nolayers");
    write_dataset(&input, 0.2);
    fs::remove_file(input.join("spliced.mtx")).unwrap();
    fs::remove_file(input.join("unspliced.mtx")).unwrap();
    let err = build_reference_file(config(&input, dir.path(), None), &dir.path().join("ref.json"))
        .unwrap_err()
        .to_string();
    assert!(err.contains("layers"), "{err}");
}
