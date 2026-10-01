//! Reader for splice-junction count matrices (input level L2).
//!
//! Supported source: a directory with `matrix.mtx[.gz]` (junction x cell),
//! `features.tsv[.gz]` describing the junctions and `barcodes.tsv[.gz]`, as
//! written by STARsolo `--soloFeatures SJ` (`Solo.out/SJ/{raw,filtered}`).
//! Feature lines are either the STARsolo / `SJ.out.tab` layout
//! `chrom  start  end  strand  motif  annotated` (strand 0 = unknown, 1 = +,
//! 2 = -; `annotated` 0/1; extra columns ignored) or a single
//! `chrom:start-end[:strand]` token, optionally followed by an `annotated`
//! column. Coordinates are the 1-based first and last intronic bases, as in
//! STAR. Cells are matched to the main matrix by barcode when the directory
//! carries `barcodes.tsv`, positionally otherwise.

use std::collections::HashMap;
use std::io::BufRead;
use std::path::{Path, PathBuf};

use tracing::{info, warn};

use crate::input::error::InputError;
use crate::io::layers::parse_mtx_counts;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
    Unknown,
}

impl Strand {
    pub fn as_char(&self) -> char {
        match self {
            Strand::Plus => '+',
            Strand::Minus => '-',
            Strand::Unknown => '.',
        }
    }
}

/// One splice junction (intron) with 1-based inclusive intron bounds.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Junction {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    /// Present in the annotation the aligner used.
    pub annotated: bool,
}

impl Junction {
    /// Donor (5' splice site) coordinate: intron start on `+`, intron end on `-`.
    pub fn donor(&self) -> u64 {
        match self.strand {
            Strand::Minus => self.end,
            _ => self.start,
        }
    }

    /// Acceptor (3' splice site) coordinate.
    pub fn acceptor(&self) -> u64 {
        match self.strand {
            Strand::Minus => self.start,
            _ => self.end,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum JunctionLocation {
    MtxDir(PathBuf),
}

impl JunctionLocation {
    pub fn describe(&self) -> String {
        match self {
            JunctionLocation::MtxDir(p) => format!("junction-mtx-dir:{}", p.display()),
        }
    }
}

#[derive(Debug, Clone)]
pub struct RawJunctions {
    pub junctions: Vec<Junction>,
    /// `(junction, cell, count)` in the main matrix's original cell index space.
    pub triplets: Vec<(u32, u32, u32)>,
    pub source: String,
    pub cells_without_junctions: usize,
}

/// Finds a junction matrix for a 10x-style input directory: an explicit
/// override, then the STARsolo sibling `<root>/SJ/<subset>` of
/// `<root>/Gene/<subset>`.
pub fn detect_junctions(
    input_dir: &Path,
    override_path: Option<&Path>,
) -> Option<JunctionLocation> {
    if let Some(p) = override_path {
        return Some(JunctionLocation::MtxDir(p.to_path_buf()));
    }
    if let (Some(subset), Some(gene_dir)) = (input_dir.file_name(), input_dir.parent())
        && gene_dir
            .file_name()
            .is_some_and(|n| n == "Gene" || n == "GeneFull")
        && let Some(root) = gene_dir.parent()
    {
        let candidate = root.join("SJ").join(subset);
        if is_junction_dir(&candidate) {
            return Some(JunctionLocation::MtxDir(candidate));
        }
    }
    None
}

fn is_junction_dir(dir: &Path) -> bool {
    dir.is_dir()
        && kira_scio::exists_plain_or_gz(&dir.join("matrix.mtx"))
        && kira_scio::exists_plain_or_gz(&dir.join("features.tsv"))
}

fn resolve(dir: &Path, base: &str) -> Option<PathBuf> {
    let plain = dir.join(base);
    if plain.is_file() {
        return Some(plain);
    }
    let gz = kira_scio::gz_path(&plain);
    if gz.is_file() { Some(gz) } else { None }
}

pub fn read_junctions(
    location: &JunctionLocation,
    barcodes: &[String],
) -> Result<RawJunctions, InputError> {
    let JunctionLocation::MtxDir(dir) = location;
    let matrix_path = resolve(dir, "matrix.mtx")
        .ok_or_else(|| InputError::MissingFile(dir.join("matrix.mtx").display().to_string()))?;
    let features_path = resolve(dir, "features.tsv")
        .ok_or_else(|| InputError::MissingFile(dir.join("features.tsv").display().to_string()))?;
    let junctions = parse_features(&features_path)?;

    let layer_barcodes = match resolve(dir, "barcodes.tsv") {
        Some(p) => Some(read_first_column(&p)?),
        None => None,
    };
    let (cell_map, cells_without) = match layer_barcodes {
        Some(ref lb) if lb != barcodes => {
            let index: HashMap<&str, u32> = barcodes
                .iter()
                .enumerate()
                .map(|(i, b)| (b.as_str(), i as u32))
                .collect();
            let map: Vec<Option<u32>> = lb.iter().map(|b| index.get(b.as_str()).copied()).collect();
            let matched = map.iter().filter(|m| m.is_some()).count();
            if matched == 0 {
                return Err(InputError::LayerMismatch(format!(
                    "no barcode of {} matches the main matrix",
                    dir.display()
                )));
            }
            let missing = barcodes.len().saturating_sub(matched);
            if missing > 0 {
                warn!(
                    junction_dir = %dir.display(),
                    cells_without_junctions = missing,
                    "some main-matrix cells have no column in the junction matrix"
                );
            }
            (Some(map), missing)
        }
        _ => (None, 0),
    };
    let n_layer_cells = cell_map.as_ref().map_or(barcodes.len(), |m| m.len());

    let (rows, cols, triplets) = parse_mtx_counts(&matrix_path)?;
    if rows != junctions.len() {
        return Err(InputError::LayerMismatch(format!(
            "{}: {} junction rows, features.tsv has {}",
            matrix_path.display(),
            rows,
            junctions.len()
        )));
    }
    if cols != n_layer_cells {
        return Err(InputError::LayerMismatch(format!(
            "{}: {} cells, expected {}",
            matrix_path.display(),
            cols,
            n_layer_cells
        )));
    }
    let triplets: Vec<(u32, u32, u32)> = match &cell_map {
        None => triplets,
        Some(map) => triplets
            .into_iter()
            .filter_map(|(j, c, k)| map[c as usize].map(|mc| (j, mc, k)))
            .collect(),
    };
    info!(
        junction_dir = %dir.display(),
        junctions = junctions.len(),
        annotated = junctions.iter().filter(|j| j.annotated).count(),
        nnz = triplets.len(),
        "junction matrix loaded"
    );
    Ok(RawJunctions {
        junctions,
        triplets,
        source: location.describe(),
        cells_without_junctions: cells_without,
    })
}

fn read_first_column(path: &Path) -> Result<Vec<String>, InputError> {
    let reader = kira_scio::open_maybe_gz_existing(path)
        .map_err(|e| InputError::UnsupportedInput(e.message))?;
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| InputError::io(path, e))?;
        let t = line.trim();
        if !t.is_empty() {
            out.push(t.split('\t').next().unwrap_or(t).to_string());
        }
    }
    Ok(out)
}

fn parse_features(path: &Path) -> Result<Vec<Junction>, InputError> {
    let reader = kira_scio::open_maybe_gz_existing(path)
        .map_err(|e| InputError::UnsupportedInput(e.message))?;
    let mut out = Vec::new();
    for (line_no, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| InputError::io(path, e))?;
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = t.split('\t').collect();
        let junction = parse_feature_line(&cols).ok_or_else(|| {
            InputError::InvalidJunctionFeatures(format!(
                "{}:{}: {t:?}",
                path.display(),
                line_no + 1
            ))
        })?;
        out.push(junction);
    }
    if out.is_empty() {
        return Err(InputError::InvalidJunctionFeatures(format!(
            "{}: no junctions",
            path.display()
        )));
    }
    Ok(out)
}

fn parse_strand(s: &str) -> Option<Strand> {
    match s {
        "1" | "+" => Some(Strand::Plus),
        "2" | "-" => Some(Strand::Minus),
        "0" | "." | "?" => Some(Strand::Unknown),
        _ => None,
    }
}

fn parse_annotated(s: &str) -> Option<bool> {
    match s.to_ascii_lowercase().as_str() {
        "1" | "true" | "annotated" | "yes" => Some(true),
        "0" | "false" | "novel" | "unannotated" | "no" => Some(false),
        _ => None,
    }
}

fn parse_feature_line(cols: &[&str]) -> Option<Junction> {
    // `chrom:start-end[:strand]` [annotated]
    if cols[0].contains(':') {
        let mut parts = cols[0].split(':');
        let chrom = parts.next()?.to_string();
        let range = parts.next()?;
        let (s, e) = range.split_once('-')?;
        let strand = parts
            .next()
            .map(parse_strand)
            .unwrap_or(Some(Strand::Unknown))?;
        let annotated = cols
            .get(1)
            .map(|a| parse_annotated(a))
            .unwrap_or(Some(false))?;
        return Some(Junction {
            chrom,
            start: s.parse().ok()?,
            end: e.parse().ok()?,
            strand,
            annotated,
        });
    }
    // STARsolo / SJ.out.tab: chrom start end strand motif annotated ...
    if cols.len() < 4 {
        return None;
    }
    let strand = parse_strand(cols[3])?;
    let annotated = match cols.len() {
        4 => false,
        5 => parse_annotated(cols[4])?,
        _ => parse_annotated(cols[5])?,
    };
    Some(Junction {
        chrom: cols[0].to_string(),
        start: cols[1].parse().ok()?,
        end: cols[2].parse().ok()?,
        strand,
        annotated,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_both_feature_layouts() {
        let j = parse_feature_line(&["chr1", "14830", "14969", "2", "2", "1"]).unwrap();
        assert_eq!(j.strand, Strand::Minus);
        assert!(j.annotated);
        assert_eq!(j.donor(), 14969);
        assert_eq!(j.acceptor(), 14830);
        let j = parse_feature_line(&["chr2:100-200:+", "0"]).unwrap();
        assert_eq!(j.strand, Strand::Plus);
        assert!(!j.annotated);
        assert_eq!(j.donor(), 100);
        assert_eq!(j.acceptor(), 200);
        let j = parse_feature_line(&["chrX:5-9"]).unwrap();
        assert_eq!(j.strand, Strand::Unknown);
        assert!(parse_feature_line(&["chr1", "x", "2", "1"]).is_none());
    }
}
