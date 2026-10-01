//! Tier-1 validation end to end: simulate -> run -> validate. Every spiked
//! effect must be recovered with AUROC >= 0.95 and its flag must keep the
//! false-positive rate at or below 1 % while recalling most positives.

use std::path::Path;

use kira_spliceqc::cli::config::{AnalysisMode, RunConfig, RunMode};
use kira_spliceqc::cli::run::run_pipeline;
use kira_spliceqc::validation::evaluate::{default_pairs, evaluate};
use kira_spliceqc::validation::simulate::{Effect, SimulationConfig, simulate};
use tempfile::tempdir;

#[test]
fn spiked_effects_are_recovered() {
    let dir = tempdir().unwrap();
    let sim = dir.path().join("sim");
    let config = SimulationConfig {
        n_cells: 1200,
        n_filler_genes: 400,
        n_junction_genes: 60,
        seed: 0x1234,
        ..SimulationConfig::default()
    };
    let effects = simulate(&config, &sim).unwrap();
    let count = |e: Effect| effects.iter().filter(|x| **x == e).count();
    assert_eq!(count(Effect::Cryptic), 60);
    assert_eq!(count(Effect::Ir), 60);
    assert_eq!(count(Effect::Damaged), 36);
    assert_eq!(count(Effect::Skip), 36);

    let out = dir.path().join("run");
    run_pipeline(RunConfig {
        input: sim.clone(),
        out_dir: out.clone(),
        cache_path: None,
        layers: None,
        junctions: Some(sim.join("sj")),
        metadata: None,
        stratify_by: None,
        reference: None,
        catalog: None,
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

    let report = evaluate(
        &out.join("kira-spliceqc").join("cells.tsv"),
        &sim.join("truth.tsv"),
        &default_pairs(),
    )
    .unwrap();
    let json = serde_json::to_string(&report).unwrap();
    assert!(json.contains("auroc"));
    let md = report.to_markdown();
    assert!(md.contains("| truth_cryptic |"));

    for r in report.results.iter().filter(|r| r.stratum == "all") {
        let auroc = r.auroc.unwrap_or_else(|| panic!("{}: no AUROC", r.truth));
        assert!(auroc >= 0.95, "{}: AUROC {auroc:.3}", r.truth);
        let fpr = r.flag_fpr.unwrap_or_else(|| panic!("{}: no flag FPR", r.truth));
        assert!(fpr <= 0.01, "{}: FPR {fpr:.3}", r.truth);
        let recall = r.flag_recall.unwrap();
        assert!(recall >= 0.9, "{}: recall {recall:.3}", r.truth);
        assert!(r.n_positive > 0 && r.n_negative > 0);
    }
    // Both strata are scored too.
    assert!(report.results.iter().any(|r| r.stratum == "TypeA"));
    assert!(report.results.iter().any(|r| r.stratum == "TypeB"));
    assert!(Path::new(&report.cells_file).ends_with("cells.tsv"));
}
