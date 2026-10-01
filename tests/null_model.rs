//! Null-model regression test.
//!
//! Generates a synthetic 10x dataset with **no splicing biology**: every cell
//! is exchangeable up to its library size, and counts are Poisson draws from
//! one shared gene-abundance vector. Any structure the tool reports on this
//! data is an artefact of the method.
//!
//! Phase 1 targets: every non-experimental flag fires on at most 1 % of cells
//! and |Spearman(metric, libsize)| <= 0.1 for every standardized metric.
//! The experimental composite flags are calibrated with the same deviation
//! rule as the production flags since Phase 2, so they share the 1 % bound;
//! the SIS classes (fixed cut-offs) keep a baseline bound. The dataset also
//! carries spliced/unspliced layers with one shared unspliced ratio, so the
//! Tier A metrics are tested on the same null.

use std::collections::HashMap;
use std::fs;
use std::path::Path;

use kira_spliceqc::cli::config::{AnalysisMode, RunConfig, RunMode};
use kira_spliceqc::cli::run::run_pipeline;
use kira_spliceqc::metrics::cell_cycle::{G2M_GENES, S_GENES};
use kira_spliceqc::metrics::splicing_instability::panels::{
    CONFLICT_RISK_PANEL, NMD_PANEL, RLOOP_RESOLUTION_PANEL, SPLICEOSOME_PANEL,
    SPLICING_RBP_PANEL,
};
use tempfile::tempdir;

const CATALOG: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/resources/genesets/splicing_genesets.tsv"
));

const N_CELLS: usize = 1500;
const N_FILLER_GENES: usize = 800;

/// Phase 1 targets for production (non-experimental) outputs.
const MAX_FLAG_FRACTION: f64 = 0.01;
const MAX_ABS_SPEARMAN_LIBSIZE: f64 = 0.10;
/// Raw control-corrected panel cores (not standardized) keep a small
/// residual depth dependence because mean-matched controls are an
/// approximate background; their standardized forms (SOS etc.) meet 0.1.
const MAX_ABS_SPEARMAN_RAW_CORE: f64 = 0.15;
/// Depth-binned deviations have a zero median in every library-size bin by
/// construction; a rank correlation across bins then only reflects shape
/// differences between bins, so they are checked per depth quintile instead:
/// |median| below this and a flag rate below `MAX_FLAG_FRACTION` in every
/// quintile.
const MAX_ABS_QUINTILE_MEDIAN: f64 = 0.15;
/// Experimental composite flags share the production bound since their
/// recalibration; the SIS classes keep a baseline (measured 2.9-4.3 %
/// Impaired+Broken, |rho(sis, libsize)| up to 0.18 from the raw entropy).
const MAX_EXPERIMENTAL_FLAG_FRACTION: f64 = 0.01;
const MAX_FAILURE_FRACTION: f64 = 0.08;
const MAX_EXPERIMENTAL_ABS_SPEARMAN: f64 = 0.30;
/// Shared unspliced ratio of every gene in the null layers.
const NULL_UNSPLICED_RATIO: f64 = 0.25;

// ---------------------------------------------------------------------------
// Deterministic PRNG (xorshift64*) and samplers; no external crates.
// ---------------------------------------------------------------------------

struct Rng(u64);

impl Rng {
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn gauss(&mut self) -> f64 {
        let u1 = self.uniform().max(1e-12);
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
    fn poisson(&mut self, lambda: f64) -> u32 {
        if lambda <= 0.0 {
            return 0;
        }
        if lambda < 30.0 {
            let l = (-lambda).exp();
            let mut k = 0u32;
            let mut p = 1.0;
            loop {
                p *= self.uniform();
                if p <= l {
                    return k;
                }
                k += 1;
            }
        }
        (lambda + lambda.sqrt() * self.gauss()).round().max(0.0) as u32
    }
}

fn panel_symbols() -> Vec<String> {
    let mut symbols: Vec<String> = Vec::new();
    for line in CATALOG.lines() {
        let line = line.trim_start_matches('\u{feff}').trim();
        if line.is_empty() || line.starts_with('#') || line.starts_with("geneset_id") {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        symbols.push(cols[2].to_string());
    }
    for panel in [
        SPLICEOSOME_PANEL,
        SPLICING_RBP_PANEL,
        RLOOP_RESOLUTION_PANEL,
        CONFLICT_RISK_PANEL,
        NMD_PANEL,
        // cell-cycle genes are ordinary null genes here (no cycling structure)
        S_GENES,
        G2M_GENES,
    ] {
        symbols.extend(panel.iter().map(|s| s.to_string()));
    }
    symbols.sort();
    symbols.dedup();
    symbols
}

/// Writes a null 10x dataset; returns per-cell true library sizes by barcode.
/// Junction matrix with no structure: `N_JUNCTION_GENES` three-exon genes,
/// every cell draws ~`JUNCTION_UMIS_PER_CELL` junction UMIs with the same
/// per-junction probabilities (canonical 0.44 each, skip 0.06, cryptic 0.06),
/// written in the STARsolo SJ layout next to the main matrix.
fn write_null_junctions(dir: &Path, rng: &mut Rng) {
    const N_JUNCTION_GENES: usize = 80;
    const JUNCTION_UMIS_PER_CELL: f64 = 600.0;
    let mut features = String::new();
    for g in 0..N_JUNCTION_GENES {
        let o = (g * 10_000) as u64;
        features.push_str(&format!("chr1\t{}\t{}\t1\t1\t1\n", o + 101, o + 200));
        features.push_str(&format!("chr1\t{}\t{}\t1\t1\t1\n", o + 301, o + 400));
        features.push_str(&format!("chr1\t{}\t{}\t1\t1\t0\n", o + 101, o + 400));
        features.push_str(&format!("chr1\t{}\t{}\t1\t1\t0\n", o + 101, o + 180));
    }
    let probs = [0.44, 0.44, 0.06, 0.06];
    let per_gene = JUNCTION_UMIS_PER_CELL / N_JUNCTION_GENES as f64;
    let mut entries: Vec<(usize, usize, u32)> = Vec::new();
    for cell in 0..N_CELLS {
        // junction depth scales with the cell's library size proxy
        let scale = (0.5 + rng.uniform()).max(0.3);
        for g in 0..N_JUNCTION_GENES {
            for (k, p) in probs.iter().enumerate() {
                let n = rng.poisson(per_gene * p * scale);
                if n > 0 {
                    entries.push((g * 4 + k + 1, cell + 1, n));
                }
            }
        }
    }
    let sj = dir.join("sj");
    fs::create_dir_all(&sj).unwrap();
    fs::write(sj.join("features.tsv"), features).unwrap();
    fs::write(
        sj.join("barcodes.tsv"),
        (0..N_CELLS).map(|c| format!("CELL{c:05}\n")).collect::<String>(),
    )
    .unwrap();
    let mut mtx = String::from("%%MatrixMarket matrix coordinate integer general\n");
    mtx.push_str(&format!("{} {} {}\n", N_JUNCTION_GENES * 4, N_CELLS, entries.len()));
    for (j, c, k) in &entries {
        mtx.push_str(&format!("{j} {c} {k}\n"));
    }
    fs::write(sj.join("matrix.mtx"), mtx).unwrap();
}

fn write_null_tenx(dir: &Path, seed: u64) -> HashMap<String, f64> {
    let mut rng = Rng(seed | 1);
    let mut genes = panel_symbols();
    let n_panel = genes.len();
    genes.extend((0..N_FILLER_GENES).map(|i| format!("FILLER{i:04}")));

    // Shared abundance vector: one biology for every cell.
    let weights: Vec<f64> = (0..genes.len())
        .map(|g| {
            if g < n_panel {
                (0.5 + 0.8 * rng.gauss()).exp()
            } else {
                (1.5 * rng.gauss()).exp()
            }
        })
        .collect();
    let total: f64 = weights.iter().sum();

    let mut triplets: Vec<(usize, usize, u32)> = Vec::new();
    let mut spliced: Vec<(usize, usize, u32)> = Vec::new();
    let mut unspliced: Vec<(usize, usize, u32)> = Vec::new();
    let mut libsizes = HashMap::with_capacity(N_CELLS);
    for cell in 0..N_CELLS {
        let libsize = (4000f64.ln() + 0.5 * rng.gauss()).exp().max(300.0).round();
        libsizes.insert(format!("CELL{cell:05}"), libsize);
        for (g, w) in weights.iter().enumerate() {
            let k = rng.poisson(libsize * w / total);
            if k > 0 {
                triplets.push((g + 1, cell + 1, k));
                // Layers: every UMI is unspliced with the same probability.
                let u = (0..k).filter(|_| rng.uniform() < NULL_UNSPLICED_RATIO).count() as u32;
                if u > 0 {
                    unspliced.push((g + 1, cell + 1, u));
                }
                if k - u > 0 {
                    spliced.push((g + 1, cell + 1, k - u));
                }
            }
        }
    }

    let write_mtx = |name: &str, entries: &[(usize, usize, u32)]| {
        let mut mtx = String::new();
        mtx.push_str("%%MatrixMarket matrix coordinate integer general\n");
        mtx.push_str(&format!("{} {} {}\n", genes.len(), N_CELLS, entries.len()));
        for (g, c, k) in entries {
            mtx.push_str(&format!("{g} {c} {k}\n"));
        }
        fs::write(dir.join(name), mtx).unwrap();
    };
    write_mtx("matrix.mtx", &triplets);
    write_mtx("spliced.mtx", &spliced);
    write_mtx("unspliced.mtx", &unspliced);
    write_null_junctions(dir, &mut rng);

    let features = genes
        .iter()
        .enumerate()
        .map(|(i, g)| format!("ENSG{i:08}\t{g}\tGene Expression\n"))
        .collect::<String>();
    fs::write(dir.join("features.tsv"), features).unwrap();

    let barcodes = (0..N_CELLS)
        .map(|c| format!("CELL{c:05}\n"))
        .collect::<String>();
    fs::write(dir.join("barcodes.tsv"), barcodes).unwrap();

    libsizes
}

// ---------------------------------------------------------------------------
// Statistics helpers.
// ---------------------------------------------------------------------------

fn ranks(values: &[f64]) -> Vec<f64> {
    let mut order: Vec<usize> = (0..values.len()).collect();
    order.sort_by(|&a, &b| values[a].partial_cmp(&values[b]).unwrap());
    let mut r = vec![0.0; values.len()];
    let mut i = 0;
    while i < order.len() {
        let mut j = i;
        while j + 1 < order.len() && values[order[j + 1]] == values[order[i]] {
            j += 1;
        }
        let avg = (i + j) as f64 / 2.0 + 1.0;
        for &idx in &order[i..=j] {
            r[idx] = avg;
        }
        i = j + 1;
    }
    r
}

fn spearman(x: &[f64], y: &[f64]) -> f64 {
    let rx = ranks(x);
    let ry = ranks(y);
    let n = rx.len() as f64;
    let mx = rx.iter().sum::<f64>() / n;
    let my = ry.iter().sum::<f64>() / n;
    let cov: f64 = rx.iter().zip(&ry).map(|(a, b)| (a - mx) * (b - my)).sum();
    let vx: f64 = rx.iter().map(|a| (a - mx).powi(2)).sum();
    let vy: f64 = ry.iter().map(|b| (b - my).powi(2)).sum();
    cov / (vx * vy).sqrt()
}

struct Table {
    header: Vec<String>,
    rows: Vec<Vec<String>>,
}

impl Table {
    fn read(path: &Path) -> Self {
        let text = fs::read_to_string(path).unwrap();
        let mut lines = text.lines();
        let header = lines
            .next()
            .unwrap()
            .split('\t')
            .map(str::to_string)
            .collect();
        let rows = lines
            .map(|l| l.split('\t').map(str::to_string).collect())
            .collect();
        Self { header, rows }
    }
    fn col(&self, name: &str) -> usize {
        self.header
            .iter()
            .position(|h| h == name)
            .unwrap_or_else(|| panic!("missing column {name}"))
    }
    fn f64_col(&self, name: &str) -> Vec<Option<f64>> {
        let c = self.col(name);
        self.rows.iter().map(|r| r[c].parse::<f64>().ok()).collect()
    }
    fn str_col(&self, name: &str) -> Vec<&str> {
        let c = self.col(name);
        self.rows.iter().map(|r| r[c].as_str()).collect()
    }
}

fn fraction(flags: &[&str], truthy: &str) -> f64 {
    flags.iter().filter(|v| **v == truthy).count() as f64 / flags.len() as f64
}

// ---------------------------------------------------------------------------
// Tests.
// ---------------------------------------------------------------------------

#[test]
fn null_model_baseline() {
    let input = tempdir().unwrap();
    let libsizes = write_null_tenx(input.path(), 0x5EED);
    let out = tempdir().unwrap();

    run_pipeline(RunConfig {
        input: input.path().to_path_buf(),
        out_dir: out.path().to_path_buf(),
        cache_path: None,
        layers: None,
        junctions: Some(input.path().join("sj")),
        metadata: None,
        stratify_by: None,
        reference: None,
        min_counts: 0,
        min_genes: 0,
        mode: AnalysisMode::Cell,
        run_mode: RunMode::Standalone,
        output_json: false,
        output_tsv: true,
        extended: true,
        threads: None,
        experimental_signatures: true,
    })
    .unwrap();

    let cells = Table::read(&out.path().join("cells.tsv"));
    assert_eq!(cells.rows.len(), N_CELLS);

    // 0. Tier B is present and mostly defined on this null (junction depth ~600).
    let junction_undefined = cells
        .f64_col("cryptic_3ss_fraction")
        .iter()
        .filter(|v| v.is_none())
        .count();
    assert!(junction_undefined * 20 <= N_CELLS, "{junction_undefined} cells without Tier B");

    // 1. No missing-value flood: the null data resolves every panel and layer.
    for metric in [
        "sis",
        "regulator_entropy_expr",
        "spliceosome_imbalance_expr",
        "missplicing_burden_expr",
        "spliceosome_core_expr",
        "SOS",
        "RLR",
        "SII",
        "unspliced_fraction",
        "intron_retention_index",
    ] {
        let undefined = cells.f64_col(metric).iter().filter(|v| v.is_none()).count();
        assert!(
            undefined * 100 <= N_CELLS,
            "{metric}: {undefined} undefined cells on null data"
        );
    }

    // 2a. Production flags on pure noise: at most 1 %.
    let mut report = String::new();
    for flag in [
        "nuclear_fraction_flag",
        "intron_retention_high",
        "cryptic_3ss_high",
        "exon_skip_high",
        "splice_site_shift_high",
    ] {
        let f = fraction(&cells.str_col(flag), "true");
        report.push_str(&format!("{flag}: {:.3}\n", f));
        assert!(f <= MAX_FLAG_FRACTION, "{flag} fraction {f:.3} on null data");
    }
    // 2b. Experimental composite flags: baseline bound until Phase 2.
    for flag in [
        "splice_overload_high",
        "rloop_risk_high",
        "splicing_instability_high",
        "genome_instability_splicing_flag",
    ] {
        let f = fraction(&cells.str_col(flag), "true");
        report.push_str(&format!("{flag} (experimental): {:.3}\n", f));
        assert!(f <= MAX_EXPERIMENTAL_FLAG_FRACTION, "{flag} fraction {f:.3} on null data");
    }
    let class = cells.str_col("class");
    let failure = fraction(&class, "Impaired") + fraction(&class, "Broken");
    report.push_str(&format!("Impaired+Broken: {failure:.3}\n"));
    assert!(
        failure <= MAX_FAILURE_FRACTION,
        "Impaired+Broken fraction {failure:.3} on null data"
    );

    // 3a. Library-size confounding of production metrics: |rho| <= 0.1.
    // (`regulator_entropy_expr` is the raw entropy and is depth-driven by
    // construction; only its depth-standardized form feeds the composites.)
    let names = cells.str_col("cell_name");
    let lib: Vec<f64> = names.iter().map(|n| libsizes[*n]).collect();
    let rho_with_libsize = |metric: &str| -> f64 {
        let values = cells.f64_col(metric);
        let (x, y): (Vec<f64>, Vec<f64>) = values
            .iter()
            .zip(&lib)
            .filter_map(|(v, l)| v.map(|v| (v, *l)))
            .unzip();
        spearman(&x, &y)
    };
    for (metric, bound) in [
        ("spliceosome_imbalance_expr", MAX_ABS_SPEARMAN_LIBSIZE),
        ("missplicing_burden_expr", MAX_ABS_SPEARMAN_LIBSIZE),
        ("SOS", MAX_ABS_SPEARMAN_LIBSIZE),
        ("unspliced_fraction", MAX_ABS_SPEARMAN_LIBSIZE),
        ("unspliced_fraction_dev", MAX_ABS_SPEARMAN_LIBSIZE),
        ("spliceosome_core_expr", MAX_ABS_SPEARMAN_RAW_CORE),
        ("s_score_expr", MAX_ABS_SPEARMAN_RAW_CORE),
    ] {
        let rho = rho_with_libsize(metric);
        report.push_str(&format!("spearman({metric}, libsize): {rho:+.3}\n"));
        assert!(
            rho.abs() <= bound,
            "{metric}: |spearman with libsize| = {:.3} on null data (bound {bound})",
            rho.abs()
        );
    }
    // Cell-cycle phase calls on null data (Seurat rule; reported, not bounded:
    // any positive noise counts as S/G2M by that rule).
    let cycling = fraction(&cells.str_col("cycling"), "true");
    report.push_str(&format!("cycling (Seurat rule, null): {cycling:.3}\n"));

    // 3b. Experimental composite: baseline bound.
    let rho = rho_with_libsize("sis");
    report.push_str(&format!("spearman(sis, libsize) (experimental): {rho:+.3}\n"));
    assert!(rho.abs() <= MAX_EXPERIMENTAL_ABS_SPEARMAN, "sis: |spearman| = {:.3}", rho.abs());

    // 3c. Depth-binned deviations: per depth quintile, |median| small and the
    // flag rate at most MAX_FLAG_FRACTION.
    let mut order: Vec<usize> = (0..lib.len()).collect();
    order.sort_by(|&a, &b| lib[a].partial_cmp(&lib[b]).unwrap());
    let quintile = order.len() / 5;
    for (metric, flag) in [
        ("intron_retention_index_dev", "intron_retention_high"),
        ("unspliced_fraction_dev", "nuclear_fraction_flag"),
    ] {
        let values = cells.f64_col(metric);
        let flags = cells.str_col(flag);
        for q in 0..5 {
            let sel = &order[q * quintile..(q + 1) * quintile];
            let mut v: Vec<f64> = sel.iter().filter_map(|&i| values[i]).collect();
            v.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let med = v[v.len() / 2];
            let rate = sel.iter().filter(|&&i| flags[i] == "true").count() as f64 / sel.len() as f64;
            report.push_str(&format!("{metric} quintile {q}: median {med:+.3}, {flag} rate {rate:.3}\n"));
            assert!(med.abs() <= MAX_ABS_QUINTILE_MEDIAN, "{metric} quintile {q}: median {med:+.3}");
            assert!(rate <= MAX_FLAG_FRACTION, "{flag} quintile {q}: rate {rate:.3}");
        }
    }
    println!("null-model report\n{report}");
}

#[test]
fn null_model_is_deterministic() {
    let input = tempdir().unwrap();
    write_null_tenx(input.path(), 0x5EED);
    let out1 = tempdir().unwrap();
    let out2 = tempdir().unwrap();
    for out in [&out1, &out2] {
        run_pipeline(RunConfig {
            input: input.path().to_path_buf(),
            out_dir: out.path().to_path_buf(),
            cache_path: None,
        layers: None,
        junctions: None,
        metadata: None,
        stratify_by: None,
        reference: None,
        min_counts: 0,
        min_genes: 0,
            mode: AnalysisMode::Cell,
            run_mode: RunMode::Standalone,
            output_json: false,
            output_tsv: true,
            extended: false,
            threads: None,
            experimental_signatures: false,
        })
        .unwrap();
    }
    let a = fs::read(out1.path().join("cells.tsv")).unwrap();
    let b = fs::read(out2.path().join("cells.tsv")).unwrap();
    assert_eq!(a, b);
}
