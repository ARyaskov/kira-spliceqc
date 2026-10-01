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
    let genes = [
        "SNRPC", "SF3A1", "SF3B1", "SFRS1", "SFRS2", "SFRS3", "HNRPA1", "ASCC3L1",
    ];
    let mut mtx = String::from("%%MatrixMarket matrix coordinate integer general\n");
    mtx.push_str(&format!("{} 2 {}\n", genes.len(), genes.len() * 2));
    for g in 1..=genes.len() {
        mtx.push_str(&format!("{g} 1 {}\n{g} 2 {}\n", g, g + 1));
    }
    fs::write(dir.join("matrix.mtx"), mtx).unwrap();
    fs::write(
        dir.join("features.tsv"),
        genes
            .iter()
            .enumerate()
            .map(|(i, g)| format!("g{i}\t{g}\n"))
            .collect::<String>(),
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
    assert_eq!(
        get("SRSF_SR").gene_ids.len(),
        3,
        "SFRS1-3 must resolve SRSF1-3"
    );
    assert!(get("SRSF_SR").missing.is_empty());
    assert_eq!(get("HNRNP").gene_ids.len(), 1);
    assert_eq!(
        get("U4_U6_CORE").gene_ids.len(),
        1,
        "ASCC3L1 must resolve SNRNP200"
    );
}

/// features.tsv with Ensembl ids and symbols the panels do not know; a user
/// catalog with an ensembl_id column resolves the genes by id.
#[test]
fn ensembl_ids_resolve_through_a_catalog_with_an_id_column() {
    use kira_spliceqc::cli::config::{AnalysisMode, RunConfig, RunMode};
    use kira_spliceqc::cli::run::run_pipeline;
    use kira_spliceqc::genesets::aliases::detect_species;

    let dir = tempdir().unwrap();
    let input = dir.path().join("data");
    fs::create_dir_all(&input).unwrap();
    // 5 panel genes under opaque symbols, 2 cells with distinct profiles.
    let ids = [
        "ENSG00000124562.9",
        "ENSG00000099995.1",
        "ENSG00000115524.17",
        "ENSG00000136450.14",
        "ENSG00000135486.3",
    ];
    let mut mtx = String::from("%%MatrixMarket matrix coordinate integer general\n5 2 10\n");
    // Non-proportional profiles (a scaled copy would give zero-MAD panels).
    for g in 1..=5 {
        mtx.push_str(&format!("{g} 1 {}\n{g} 2 {}\n", 3 + g, (g * 3) % 7 + 2));
    }
    fs::write(input.join("matrix.mtx"), mtx).unwrap();
    fs::write(
        input.join("features.tsv"),
        ids.iter()
            .enumerate()
            .map(|(i, id)| format!("{id}\tLOC{i}\tGene Expression\n"))
            .collect::<String>(),
    )
    .unwrap();
    fs::write(input.join("barcodes.tsv"), "c1\nc2\n").unwrap();

    let catalog = dir.path().join("catalog.tsv");
    fs::write(
        &catalog,
        "geneset_id\taxis\tgene_symbol\tensembl_id\n\
         U1_CORE\tCORE_SPLICEOSOME\tSNRPC\tENSG00000124562\n\
         U2_CORE\tCORE_SPLICEOSOME\tSF3A1\tENSG00000099995\n\
         U2_CORE\tCORE_SPLICEOSOME\tSF3B1\tENSG00000115524\n\
         SF3B_AXIS\tCORE_SPLICEOSOME\tSF3B1\tENSG00000115524\n\
         SRSF_SR\tREGULATORS\tSRSF1\tENSG00000136450\n\
         HNRNP\tREGULATORS\tHNRNPA1\tENSG00000135486\n\
         MINOR_U12\tMINOR_SPLICEOSOME\tZRSR2\t\n\
         NMD_SURVEILLANCE\tSURVEILLANCE\tUPF1\t\n",
    )
    .unwrap();

    // Symbols alone resolve nothing; the pipeline fails at the core gate.
    let stage0 = run_stage0(&input, RunMode::Standalone, None).unwrap();
    let out = tempdir().unwrap();
    let matrix = run_stage1(&stage0, out.path()).unwrap();
    assert_eq!(
        detect_species(&matrix),
        "human",
        "species from ENSG prefixes"
    );
    let by_symbol = load_catalog(Path::new("does-not-exist.tsv"), &matrix).unwrap();
    assert!(
        by_symbol
            .genesets
            .iter()
            .find(|g| g.id == "SRSF_SR")
            .unwrap()
            .gene_ids
            .is_empty()
    );

    // With the id column every panel gene resolves (version suffixes ignored).
    let by_id = load_catalog(&catalog, &matrix).unwrap();
    let get = |id: &str| by_id.genesets.iter().find(|g| g.id == id).unwrap();
    assert_eq!(get("SRSF_SR").gene_ids.len(), 1);
    assert_eq!(get("U2_CORE").gene_ids.len(), 2);
    assert_eq!(get("HNRNP").missing.len(), 0);
    assert_eq!(get("MINOR_U12").missing, vec!["ZRSR2"]);

    // --catalog goes through the whole run.
    let out2 = tempdir().unwrap();
    run_pipeline(RunConfig {
        input: input.clone(),
        out_dir: out2.path().to_path_buf(),
        cache_path: None,
        layers: None,
        junctions: None,
        metadata: None,
        stratify_by: None,
        reference: None,
        catalog: Some(catalog.clone()),
        min_counts: 0,
        min_genes: 0,
        mode: AnalysisMode::Cell,
        run_mode: RunMode::Pipeline,
        output_json: false,
        output_tsv: true,
        extended: false,
        threads: None,
        experimental_signatures: false,
    })
    .unwrap();
    let summary: serde_json::Value = serde_json::from_slice(
        &fs::read(out2.path().join("kira-spliceqc").join("summary.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(summary["input"]["species"], "human");
    assert!(
        summary["provenance"]["geneset_catalog"]["source"]
            .as_str()
            .unwrap()
            .ends_with("catalog.tsv")
    );
}
