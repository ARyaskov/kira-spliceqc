//! Synthetic datasets with known truth for tier-1 validation.
//!
//! Every cell is a Poisson draw from one shared gene-abundance vector (two
//! cell types differ only in their baseline unspliced fraction), so the
//! dataset carries no splicing structure except the spiked effects:
//!
//! - `cryptic`: a fraction of the cell's reads at every donor use the cryptic
//!   acceptor 20 nt upstream of the canonical one (SF3B1-like);
//! - `ir`: the unspliced counts of every gene are multiplied (global intron
//!   retention);
//! - `damaged`: the unspliced counts are divided by 10 (cytoplasmic debris);
//! - `skip`: a fraction of the reads at every three-exon gene skip the
//!   middle exon.
//!
//! A cell carries at most one effect. The truth table lists every effect as a
//! boolean column so `kira-spliceqc validate` can score metrics and flags.

use std::fs;
use std::path::Path;

use crate::input::error::InputError;
use crate::metrics::cell_cycle::{G2M_GENES, S_GENES};
use crate::metrics::splicing_instability::panels::{
    CONFLICT_RISK_PANEL, NMD_PANEL, RLOOP_RESOLUTION_PANEL, SPLICEOSOME_PANEL, SPLICING_RBP_PANEL,
};

#[derive(Debug, Clone)]
pub struct SimulationConfig {
    pub n_cells: usize,
    pub n_filler_genes: usize,
    pub n_junction_genes: usize,
    pub seed: u64,
    pub cryptic_fraction: f64,
    pub cryptic_ratio: f64,
    pub ir_fraction: f64,
    pub ir_fold: f64,
    pub damaged_fraction: f64,
    pub skip_fraction: f64,
    pub skip_ratio: f64,
}

impl Default for SimulationConfig {
    fn default() -> Self {
        Self {
            n_cells: 2000,
            n_filler_genes: 800,
            n_junction_genes: 80,
            seed: 0x5EED,
            cryptic_fraction: 0.05,
            cryptic_ratio: 0.15,
            ir_fraction: 0.05,
            ir_fold: 2.0,
            damaged_fraction: 0.03,
            skip_fraction: 0.03,
            skip_ratio: 0.25,
        }
    }
}

pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Self(seed | 1)
    }
    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    pub fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    pub fn gauss(&mut self) -> f64 {
        let u1 = self.uniform().max(1e-12);
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
    pub fn poisson(&mut self, lambda: f64) -> u32 {
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
    pub fn binomial(&mut self, n: u32, p: f64) -> u32 {
        (0..n).filter(|_| self.uniform() < p).count() as u32
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Effect {
    None,
    Cryptic,
    Ir,
    Damaged,
    Skip,
}

impl Effect {
    fn column(&self) -> Option<&'static str> {
        match self {
            Effect::None => None,
            Effect::Cryptic => Some("truth_cryptic"),
            Effect::Ir => Some("truth_ir"),
            Effect::Damaged => Some("truth_damaged"),
            Effect::Skip => Some("truth_skip"),
        }
    }
}

/// Splicing-panel genes that must exist for the expression stages to run.
pub fn panel_symbols() -> Vec<String> {
    const CATALOG: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/resources/genesets/splicing_genesets.tsv"));
    let mut symbols: Vec<String> = Vec::new();
    for line in CATALOG.lines() {
        let line = line.trim_start_matches('\u{feff}').trim();
        if line.is_empty() || line.starts_with('#') || line.starts_with("geneset_id") {
            continue;
        }
        if let Some(sym) = line.split('\t').nth(2) {
            symbols.push(sym.to_string());
        }
    }
    for panel in [SPLICEOSOME_PANEL, SPLICING_RBP_PANEL, RLOOP_RESOLUTION_PANEL, CONFLICT_RISK_PANEL, NMD_PANEL, S_GENES, G2M_GENES] {
        symbols.extend(panel.iter().map(|s| s.to_string()));
    }
    symbols.sort();
    symbols.dedup();
    symbols
}

/// Writes the dataset (10x directory with layers, `sj/` junction matrix,
/// `metadata.tsv`, `truth.tsv`) and returns the effect of every cell.
pub fn simulate(config: &SimulationConfig, out: &Path) -> Result<Vec<Effect>, InputError> {
    fs::create_dir_all(out).map_err(|e| InputError::io(out, e))?;
    let mut rng = Rng::new(config.seed);
    let n_cells = config.n_cells;

    // Genes and shared abundances.
    let mut genes = panel_symbols();
    let n_panel = genes.len();
    genes.extend((0..config.n_filler_genes).map(|i| format!("FILLER{i:04}")));
    let weights: Vec<f64> = (0..genes.len())
        .map(|g| if g < n_panel { (0.5 + 0.8 * rng.gauss()).exp() } else { (1.5 * rng.gauss()).exp() })
        .collect();
    let total: f64 = weights.iter().sum();

    // Effects: at most one per cell, assigned deterministically from the seed.
    let mut effects = vec![Effect::None; n_cells];
    let assign = |effect: Effect, fraction: f64, effects: &mut Vec<Effect>, rng: &mut Rng| {
        let target = (fraction * n_cells as f64).round() as usize;
        let mut placed = 0;
        let mut guard = 0;
        while placed < target && guard < n_cells * 10 {
            let c = (rng.next_u64() % n_cells as u64) as usize;
            if effects[c] == Effect::None {
                effects[c] = effect;
                placed += 1;
            }
            guard += 1;
        }
    };
    assign(Effect::Cryptic, config.cryptic_fraction, &mut effects, &mut rng);
    assign(Effect::Ir, config.ir_fraction, &mut effects, &mut rng);
    assign(Effect::Damaged, config.damaged_fraction, &mut effects, &mut rng);
    assign(Effect::Skip, config.skip_fraction, &mut effects, &mut rng);

    // Two cell types with different baseline unspliced fractions.
    let cell_type: Vec<&str> = (0..n_cells).map(|c| if c % 2 == 0 { "TypeA" } else { "TypeB" }).collect();
    let base_uf = |c: usize| if c.is_multiple_of(2) { 0.25 } else { 0.45 };

    let mut matrix = Vec::new();
    let mut spliced = Vec::new();
    let mut unspliced = Vec::new();
    let mut libsizes = Vec::with_capacity(n_cells);
    for (c, effect) in effects.iter().enumerate() {
        let libsize = (4000f64.ln() + 0.5 * rng.gauss()).exp().max(300.0).round();
        libsizes.push(libsize);
        let uf = match effect {
            Effect::Ir => (base_uf(c) * config.ir_fold).min(0.95),
            Effect::Damaged => base_uf(c) * 0.1,
            _ => base_uf(c),
        };
        for (g, w) in weights.iter().enumerate() {
            let k = rng.poisson(libsize * w / total);
            if k > 0 {
                matrix.push((g + 1, c + 1, k));
                let u = rng.binomial(k, uf);
                if u > 0 {
                    unspliced.push((g + 1, c + 1, u));
                }
                if k - u > 0 {
                    spliced.push((g + 1, c + 1, k - u));
                }
            }
        }
    }
    let write_mtx = |path: &Path, n_rows: usize, entries: &[(usize, usize, u32)]| -> Result<(), InputError> {
        let mut s = String::from("%%MatrixMarket matrix coordinate integer general\n");
        s.push_str(&format!("{n_rows} {n_cells} {}\n", entries.len()));
        for (r, c, k) in entries {
            s.push_str(&format!("{r} {c} {k}\n"));
        }
        fs::write(path, s).map_err(|e| InputError::io(path, e))
    };
    write_mtx(&out.join("matrix.mtx"), genes.len(), &matrix)?;
    write_mtx(&out.join("spliced.mtx"), genes.len(), &spliced)?;
    write_mtx(&out.join("unspliced.mtx"), genes.len(), &unspliced)?;
    let features: String = genes.iter().enumerate().map(|(i, g)| format!("ENSG{i:08}\t{g}\tGene Expression\n")).collect();
    fs::write(out.join("features.tsv"), features).map_err(|e| InputError::io(out, e))?;
    let barcodes: String = (0..n_cells).map(|c| format!("CELL{c:05}\n")).collect();
    fs::write(out.join("barcodes.tsv"), &barcodes).map_err(|e| InputError::io(out, e))?;

    // Junctions: three-exon genes; cryptic acceptor and skip junctions.
    let sj = out.join("sj");
    fs::create_dir_all(&sj).map_err(|e| InputError::io(&sj, e))?;
    let mut features = String::new();
    for g in 0..config.n_junction_genes {
        let o = (g * 10_000) as u64;
        features.push_str(&format!("chr1\t{}\t{}\t1\t1\t1\n", o + 101, o + 200));
        features.push_str(&format!("chr1\t{}\t{}\t1\t1\t1\n", o + 301, o + 400));
        features.push_str(&format!("chr1\t{}\t{}\t1\t1\t0\n", o + 101, o + 400));
        features.push_str(&format!("chr1\t{}\t{}\t1\t1\t0\n", o + 101, o + 180));
    }
    fs::write(sj.join("features.tsv"), features).map_err(|e| InputError::io(&sj, e))?;
    fs::write(sj.join("barcodes.tsv"), &barcodes).map_err(|e| InputError::io(&sj, e))?;
    let per_gene = 600.0 / config.n_junction_genes as f64;
    let mut entries = Vec::new();
    for c in 0..n_cells {
        let scale = (libsizes[c] / 4000.0).clamp(0.3, 3.0);
        let (p_cryptic, p_skip) = match effects[c] {
            Effect::Cryptic => (config.cryptic_ratio, 0.03),
            Effect::Skip => (0.03, config.skip_ratio),
            _ => (0.03, 0.03),
        };
        for g in 0..config.n_junction_genes {
            let donor_reads = rng.poisson(per_gene * 0.5 * scale);
            let cryptic = rng.binomial(donor_reads, p_cryptic);
            let skip = rng.binomial(donor_reads - cryptic, p_skip);
            let canonical = donor_reads - cryptic - skip;
            let second = rng.poisson(per_gene * 0.5 * scale);
            for (k, count) in [(0usize, canonical), (1, second), (2, skip), (3, cryptic)] {
                if count > 0 {
                    entries.push((g * 4 + k + 1, c + 1, count));
                }
            }
        }
    }
    write_mtx(&sj.join("matrix.mtx"), config.n_junction_genes * 4, &entries)?;

    // Metadata (strata) and truth.
    let mut metadata = String::from("barcode\tcell_type\n");
    let mut truth = String::from("barcode\tcell_type\tlibsize_true\teffect\ttruth_cryptic\ttruth_ir\ttruth_damaged\ttruth_skip\n");
    for c in 0..n_cells {
        metadata.push_str(&format!("CELL{c:05}\t{}\n", cell_type[c]));
        let col = effects[c].column();
        let flag = |name: &str| if col == Some(name) { "true" } else { "false" };
        truth.push_str(&format!(
            "CELL{c:05}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            cell_type[c],
            libsizes[c],
            col.map(|s| s.trim_start_matches("truth_")).unwrap_or("none"),
            flag("truth_cryptic"),
            flag("truth_ir"),
            flag("truth_damaged"),
            flag("truth_skip"),
        ));
    }
    fs::write(out.join("metadata.tsv"), metadata).map_err(|e| InputError::io(out, e))?;
    fs::write(out.join("truth.tsv"), truth).map_err(|e| InputError::io(out, e))?;
    Ok(effects)
}
