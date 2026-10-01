//! End-to-end Tier B: a STARsolo layout (`Solo.out/Gene/raw` + `Solo.out/SJ/raw`)
//! with 60 genes of three exons each. Normal cells use the annotated
//! junctions; 6 "SF3B1-like" cells put 15 % of their donor reads on cryptic
//! acceptors 20 nt upstream of the canonical ones; 4 cells skip the middle
//! exon in a quarter of their reads.

use std::fs;
use std::path::Path;

use kira_spliceqc::cli::config::{AnalysisMode, RunConfig, RunMode};
use kira_spliceqc::cli::run::run_pipeline;
use kira_spliceqc::io::junctions::JunctionLocation;
use kira_spliceqc::pipeline::stage0_input::run_stage0_full;
use tempfile::tempdir;

const N_CELLS: usize = 150;
const N_GENES: usize = 60;

fn write_solo(root: &Path) {
    // --- Gene/raw: the main matrix (5 panel genes + fillers) -----------------
    let gene_dir = root.join("Gene").join("raw");
    fs::create_dir_all(&gene_dir).unwrap();
    let mut genes: Vec<String> = ["SNRPC", "SF3A1", "SF3B1", "SRSF1", "HNRNPA1"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    genes.extend((5..N_GENES).map(|i| format!("GENE{i:03}")));
    let mut m = String::new();
    let mut entries = 0;
    for c in 0..N_CELLS {
        for g in 1..=genes.len() {
            m.push_str(&format!("{g} {} {}\n", c + 1, 20 + ((c * 7 + g * 3) % 11)));
            entries += 1;
        }
    }
    fs::write(
        gene_dir.join("matrix.mtx"),
        format!(
            "%%MatrixMarket matrix coordinate integer general\n{} {} {}\n{m}",
            genes.len(),
            N_CELLS,
            entries
        ),
    )
    .unwrap();
    fs::write(
        gene_dir.join("features.tsv"),
        genes
            .iter()
            .enumerate()
            .map(|(i, g)| format!("g{i}\t{g}\n"))
            .collect::<String>(),
    )
    .unwrap();
    let barcodes: String = (0..N_CELLS).map(|c| format!("CELL{c:03}\n")).collect();
    fs::write(gene_dir.join("barcodes.tsv"), &barcodes).unwrap();

    // --- SJ/raw: 4 junctions per gene ------------------------------------------
    // gene g occupies chr1 at offset g*10000: exon1 1..100, exon2 201..300, exon3 401..500
    // j0 = 101-200 annotated, j1 = 301-400 annotated, j2 = 101-400 skip (novel),
    // j3 = 101-180 cryptic acceptor 20 nt upstream (novel).
    let sj_dir = root.join("SJ").join("raw");
    fs::create_dir_all(&sj_dir).unwrap();
    let mut features = String::new();
    for g in 0..N_GENES {
        let o = (g * 10_000) as u64;
        features.push_str(&format!("chr1\t{}\t{}\t1\t1\t1\n", o + 101, o + 200));
        features.push_str(&format!("chr1\t{}\t{}\t1\t1\t1\n", o + 301, o + 400));
        features.push_str(&format!("chr1\t{}\t{}\t1\t1\t0\n", o + 101, o + 400));
        features.push_str(&format!("chr1\t{}\t{}\t1\t1\t0\n", o + 101, o + 180));
    }
    fs::write(sj_dir.join("features.tsv"), features).unwrap();
    fs::write(sj_dir.join("barcodes.tsv"), &barcodes).unwrap();
    let mut sj = String::new();
    let mut nnz = 0;
    for c in 0..N_CELLS {
        let mutant = c < 6;
        let skipper = (6..10).contains(&c);
        for g in 0..N_GENES {
            let base = 10 + ((c * 5 + g * 3) % 7) as u32; // reads from donor 1 of this gene
            let cryptic = if mutant {
                (base as f64 * 0.15).round() as u32
            } else {
                (g + c) as u32 % 2
            };
            let skip = if skipper {
                (base as f64 * 0.25).round() as u32
            } else {
                (g * 3 + c) as u32 % 2
            };
            let j0 = base - cryptic.min(base / 2);
            let j1 = base + ((c + g) % 3) as u32;
            let rows = [(0, j0), (1, j1), (2, skip), (3, cryptic)];
            for (k, count) in rows {
                if count > 0 {
                    sj.push_str(&format!("{} {} {count}\n", g * 4 + k + 1, c + 1));
                    nnz += 1;
                }
            }
        }
    }
    fs::write(
        sj_dir.join("matrix.mtx"),
        format!(
            "%%MatrixMarket matrix coordinate integer general\n{} {} {}\n{sj}",
            N_GENES * 4,
            N_CELLS,
            nnz
        ),
    )
    .unwrap();
}

fn config(input: &Path, out: &Path) -> RunConfig {
    RunConfig {
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
        output_json: true,
        output_tsv: true,
        extended: false,
        threads: None,
        experimental_signatures: false,
    }
}

#[test]
fn starsolo_sj_sibling_is_detected() {
    let dir = tempdir().unwrap();
    let root = dir.path().join("Solo.out");
    write_solo(&root);
    let gene_dir = root.join("Gene").join("raw");
    let stage0 = run_stage0_full(&gene_dir, RunMode::Standalone, None, None, None).unwrap();
    assert_eq!(
        stage0.junctions,
        Some(JunctionLocation::MtxDir(root.join("SJ").join("raw")))
    );
    // Explicit override wins; a missing override path is an error.
    let elsewhere = root.join("SJ").join("raw");
    let stage0 =
        run_stage0_full(&gene_dir, RunMode::Standalone, None, None, Some(&elsewhere)).unwrap();
    assert_eq!(stage0.junctions, Some(JunctionLocation::MtxDir(elsewhere)));
    assert!(
        run_stage0_full(
            &gene_dir,
            RunMode::Standalone,
            None,
            None,
            Some(&root.join("nope"))
        )
        .is_err()
    );
}

#[test]
fn cryptic_and_skip_cells_are_flagged() {
    let dir = tempdir().unwrap();
    let root = dir.path().join("Solo.out");
    write_solo(&root);
    let out = tempdir().unwrap();
    run_pipeline(config(&root.join("Gene").join("raw"), out.path())).unwrap();
    let base = out.path().join("kira-spliceqc");

    let summary: serde_json::Value =
        serde_json::from_slice(&fs::read(base.join("summary.json")).unwrap()).unwrap();
    assert_eq!(summary["input"]["levels"], serde_json::json!(["L0", "L2"]));
    let j = &summary["junctions"];
    assert_eq!(j["n_junctions"], N_GENES * 4);
    assert_eq!(j["n_annotated"], N_GENES * 2);
    assert_eq!(j["n_cryptic_acceptor_junctions"], N_GENES);
    assert_eq!(j["n_skip_junctions"], N_GENES);
    assert_eq!(j["n_defined_cells"], N_CELLS);
    assert_eq!(j["cells_without_junctions"], 0);
    assert_eq!(
        summary["provenance"]["input_levels"],
        serde_json::json!(["L0", "L2"])
    );
    assert_eq!(
        summary["provenance"]["parameters"]["cryptic_window_nt"],
        serde_json::json!([10, 50])
    );

    let cells = fs::read_to_string(base.join("cells.tsv")).unwrap();
    let mut lines = cells.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let col = |n: &str| header.iter().position(|h| *h == n).unwrap();
    let (name, cf, ch, sf, sh, shift) = (
        col("cell_name"),
        col("cryptic_3ss_fraction"),
        col("cryptic_3ss_high"),
        col("exon_skip_fraction"),
        col("exon_skip_high"),
        col("splice_site_shift"),
    );
    let mut cryptic_flagged = Vec::new();
    let mut skip_flagged = Vec::new();
    for line in lines {
        let f: Vec<&str> = line.split('\t').collect();
        let idx: usize = f[name][4..].parse().unwrap();
        let cryptic: f64 = f[cf].parse().unwrap();
        let skip: f64 = f[sf].parse().unwrap();
        if idx < 6 {
            assert!(
                cryptic > 0.10,
                "mutant cell {idx}: cryptic fraction {cryptic}"
            );
        } else {
            assert!(
                cryptic < 0.08,
                "normal cell {idx}: cryptic fraction {cryptic}"
            );
        }
        if (6..10).contains(&idx) {
            assert!(skip > 0.15, "skipper cell {idx}: skip fraction {skip}");
        }
        if f[ch] == "true" {
            cryptic_flagged.push(idx);
        }
        if f[sh] == "true" {
            skip_flagged.push(idx);
        }
        assert!(!f[shift].is_empty(), "shift defined for cell {idx}");
    }
    assert_eq!(cryptic_flagged, vec![0, 1, 2, 3, 4, 5]);
    assert_eq!(skip_flagged, vec![6, 7, 8, 9]);

    let cj: serde_json::Value =
        serde_json::from_slice(&fs::read(base.join("cells.json")).unwrap()).unwrap();
    assert_eq!(cj["input_levels"], serde_json::json!(["L0", "L2"]));
    assert!(cj["cells"][0]["junctions"]["junction_umis"].is_number());
    assert_eq!(
        cj["junctions"]["n_site_groups"].as_u64().unwrap(),
        (N_GENES * 2) as u64
    );
}
