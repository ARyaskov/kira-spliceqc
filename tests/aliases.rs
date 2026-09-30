//! Legacy gene symbols resolve to the current panel genes.

use std::fs;
use std::path::Path;

use kira_spliceqc::cli::config::RunMode;
use kira_spliceqc::genesets::load_catalog;
use kira_spliceqc::pipeline::stage0_input::run_stage0;
use kira_spliceqc::pipeline::stage1_expression::run_stage1;
use tempfile::tempdir;

fn write_legacy_tenx(dir: &Path) {
    fs::create_dir_all(dir).unwrap();
    // Old annotation: SFRS1/SFRS2/SFRS3 for SRSF1-3, HNRPA1 for HNRNPA1, ASCC3L1 for SNRNP200.
    let genes = ["SNRPC", "SF3A1", "SF3B1", "SFRS1", "SFRS2", "SFRS3", "HNRPA1", "ASCC3L1"];
    let mut mtx = String::from("%%MatrixMarket matrix coordinate integer general\n");
    mtx.push_str(&format!("{} 2 {}\n", genes.len(), genes.len() * 2));
    for g in 1..=genes.len() {
        mtx.push_str(&format!("{g} 1 {}\n{g} 2 {}\n", g, g + 1));
    }
    fs::write(dir.join("matrix.mtx"), mtx).unwrap();
    fs::write(
        dir.join("features.tsv"),
        genes.iter().enumerate().map(|(i, g)| format!("g{i}\t{g}\n")).collect::<String>(),
    )
    .unwrap();
    fs::write(dir.join("barcodes.tsv"), "c1\nc2\n").unwrap();
}

#[test]
fn legacy_symbols_resolve_panel_genes() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("data");
    write_legacy_tenx(&input);
    let stage0 = run_stage0(&input, RunMode::Standalone, None).unwrap();
    let out = tempdir().unwrap();
    let matrix = run_stage1(&stage0, out.path()).unwrap();
    let catalog = load_catalog(Path::new("does-not-exist.tsv"), &matrix).unwrap();
    let get = |id: &str| catalog.genesets.iter().find(|g| g.id == id).unwrap();
    assert_eq!(get("SRSF_SR").gene_ids.len(), 3, "SFRS1-3 must resolve SRSF1-3");
    assert!(get("SRSF_SR").missing.is_empty());
    assert_eq!(get("HNRNP").gene_ids.len(), 1);
    assert_eq!(get("U4_U6_CORE").gene_ids.len(), 1, "ASCC3L1 must resolve SNRNP200");
}
