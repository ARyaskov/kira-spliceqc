use std::fs;

use kira_spliceqc::model::assembly_phase::AssemblyPhaseImbalanceMetrics;
use kira_spliceqc::model::coupling::CouplingStressMetrics;
use kira_spliceqc::model::exon_intron_bias::ExonIntronDefinitionMetrics;
use kira_spliceqc::model::imbalance::SpliceosomeImbalanceMetrics;
use kira_spliceqc::model::isoform_dispersion::IsoformDispersionMetrics;
use kira_spliceqc::model::missplicing::MissplicingMetrics;
use kira_spliceqc::model::sis::{SpliceIntegrityClass, SpliceIntegrityMetrics};
use kira_spliceqc::output::{json, summary, tsv};
use kira_spliceqc::pipeline::stage7_output::{OutputOptions, run_stage7};
use tempfile::tempdir;

fn make_metrics(
    n: usize,
) -> (
    Vec<String>,
    IsoformDispersionMetrics,
    MissplicingMetrics,
    SpliceosomeImbalanceMetrics,
    SpliceIntegrityMetrics,
    CouplingStressMetrics,
    ExonIntronDefinitionMetrics,
    AssemblyPhaseImbalanceMetrics,
) {
    let cells = (0..n).map(|i| format!("cell{}", i)).collect::<Vec<_>>();
    let isoform = IsoformDispersionMetrics {
        entropy: vec![0.9; n],
        dispersion: vec![0.8; n],
        z_entropy: vec![1.0; n],
    };
    let missplicing = MissplicingMetrics {
        b_core: vec![0.1; n],
        b_u12: vec![0.2; n],
        b_nmd: vec![0.3; n],
        b_srhn: vec![0.4; n],
        burden: vec![0.5; n],
        burden_star: vec![0.6; n],
    };
    let imbalance = SpliceosomeImbalanceMetrics {
        z_u1: vec![0.1; n],
        z_u2: vec![0.2; n],
        z_sf3b: vec![0.3; n],
        z_srsf: vec![0.4; n],
        z_hnrnp: vec![0.5; n],
        z_u12: vec![0.6; n],
        z_nmd: vec![0.7; n],
        axis_sr_hnrnp: vec![0.8; n],
        axis_u2_u1: vec![0.9; n],
        axis_u12_major: vec![1.0; n],
        axis_nmd: vec![1.1; n],
        imbalance: vec![1.2; n],
    };
    let sis = SpliceIntegrityMetrics {
        sis: vec![0.75; n],
        class: vec![SpliceIntegrityClass::Stressed; n],
        p_missplicing: vec![0.1; n],
        p_imbalance: vec![0.2; n],
        p_entropy_z: vec![0.3; n],
        p_entropy_abs: vec![0.4; n],
    };
    let coupling = CouplingStressMetrics {
        coupling_stress: vec![0.05; n],
    };
    let exon_intron = ExonIntronDefinitionMetrics {
        exon_definition_bias: vec![0.2; n],
        z_srsf: vec![0.1; n],
        z_u2af: vec![0.2; n],
        z_hnrnp: vec![0.3; n],
    };
    let assembly = AssemblyPhaseImbalanceMetrics {
        z_ea: vec![0.1; n],
        z_b: vec![0.2; n],
        z_cat: vec![0.3; n],
        ea_imbalance: vec![0.01; n],
        b_imbalance: vec![0.02; n],
        cat_imbalance: vec![0.03; n],
    };
    (
        cells,
        isoform,
        missplicing,
        imbalance,
        sis,
        coupling,
        exon_intron,
        assembly,
    )
}

#[test]
fn json_schema_sanity() {
    let (cells, isoform, missplicing, imbalance, sis, coupling, exon_intron, assembly) =
        make_metrics(2);
    let dir = tempdir().unwrap();
    let path = dir.path().join("spliceqc.json");
    json::write_json(
        &path,
        &cells,
        &isoform,
        &missplicing,
        &imbalance,
        &sis,
        Some(&coupling),
        Some(&exon_intron),
        Some(&assembly),
        None,
        None,
        None,
        None,
    )
    .unwrap();

    let data = fs::read_to_string(&path).unwrap();
    let v: serde_json::Value = serde_json::from_str(&data).unwrap();
    assert_eq!(v["schema_version"], "1.0");
    assert_eq!(v["tool"], "kira-spliceqc");
    assert_eq!(v["mode"], "cell");
    assert_eq!(v["n_cells"], 2);
    assert!(v["cells"].is_array());
    assert!(v["cells"][0]["coupling"]["coupling_stress"].is_number());
    assert!(v["cells"][0]["exon_intron_bias"]["exon_definition_bias"].is_number());
    assert!(v["cells"][0]["assembly_phase"]["ea_imbalance"].is_number());
}

#[test]
fn tsv_header_order() {
    let (cells, isoform, missplicing, imbalance, sis, coupling, exon_intron, assembly) =
        make_metrics(1);
    let dir = tempdir().unwrap();
    let path = dir.path().join("spliceqc.tsv");
    tsv::write_tsv(
        &path,
        &cells,
        &isoform,
        &missplicing,
        &imbalance,
        &sis,
        &coupling,
        &exon_intron,
        &assembly,
    )
    .unwrap();

    let data = fs::read_to_string(&path).unwrap();
    let header = data.lines().next().unwrap();
    assert_eq!(header, tsv::header());
}

#[test]
fn summary_formatting_snapshot() {
    let (_cells, _isoform, _missplicing, _imbalance, sis, _coupling, _exon_intron, _assembly) =
        make_metrics(3);
    let text = summary::format_summary(&sis, None, None, None);
    let expected = "kira-spliceqc summary\n---------------------\nCells analyzed: 3\n\nIntegrity classes:\n  Intact:       0 (0.0%)\n  Stressed:     3 (100.0%)\n  Impaired:     0 (0.0%)\n  Broken:       0 (0.0%)\n\nMedian SIS: 0.75\nFailure fraction (Impaired+Broken): 0.0%\n\nCryptic splicing risk > 0.7: N/A\nSpliceosome collapse: N/A\nCell-cycle confounded: N/A\n";
    assert_eq!(text, expected);
}

#[test]
fn json_deterministic_bytes() {
    let (cells, isoform, missplicing, imbalance, sis, coupling, exon_intron, assembly) =
        make_metrics(2);
    let dir = tempdir().unwrap();
    let path1 = dir.path().join("spliceqc1.json");
    let path2 = dir.path().join("spliceqc2.json");

    json::write_json(
        &path1,
        &cells,
        &isoform,
        &missplicing,
        &imbalance,
        &sis,
        Some(&coupling),
        Some(&exon_intron),
        Some(&assembly),
        None,
        None,
        None,
        None,
    )
    .unwrap();
    json::write_json(
        &path2,
        &cells,
        &isoform,
        &missplicing,
        &imbalance,
        &sis,
        Some(&coupling),
        Some(&exon_intron),
        Some(&assembly),
        None,
        None,
        None,
        None,
    )
    .unwrap();

    let b1 = fs::read(&path1).unwrap();
    let b2 = fs::read(&path2).unwrap();
    assert_eq!(b1, b2);
}

#[test]
fn run_stage7_outputs() {
    let (cells, isoform, missplicing, imbalance, sis, coupling, exon_intron, assembly) =
        make_metrics(1);
    let dir = tempdir().unwrap();
    let summary = run_stage7(
        dir.path(),
        &cells,
        &isoform,
        &missplicing,
        &imbalance,
        &sis,
        Some(&coupling),
        Some(&exon_intron),
        Some(&assembly),
        None,
        None,
        None,
        None,
        OutputOptions {
            json: false,
            tsv: false,
        },
    )
    .unwrap();
    assert!(summary.contains("kira-spliceqc summary"));
    assert!(dir.path().join("spliceqc.json").exists());
    assert!(dir.path().join("spliceqc.tsv").exists());
}
