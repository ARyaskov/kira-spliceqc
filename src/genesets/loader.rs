use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use ahash::AHashMap;
use tracing::{debug, info, warn};

use crate::expression::ExpressionMatrix;
use crate::genesets::{Geneset, GenesetCatalog};
use crate::input::error::InputError;

const EMBEDDED_SPLICE_GENESETS: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/resources/genesets/splicing_genesets.tsv"
));

pub fn load_catalog(
    path: &Path,
    matrix: &dyn ExpressionMatrix,
) -> Result<GenesetCatalog, InputError> {
    let mut entries: BTreeMap<String, (String, Vec<String>)> = BTreeMap::new();
    match File::open(path) {
        Ok(file) => {
            let reader = BufReader::new(file);
            parse_catalog_lines(reader.lines(), path.display().to_string(), &mut entries)?
        }
        Err(_) => parse_catalog_str(
            EMBEDDED_SPLICE_GENESETS,
            "embedded://splicing_genesets.tsv",
            &mut entries,
        )?,
    }

    // Case-insensitive index: matrix gene symbols may be mixed-case (mouse "Snrpb").
    // First occurrence wins (matches stage15 splicing_instability behavior).
    let mut symbol_to_id: AHashMap<String, u32> = AHashMap::with_capacity(matrix.n_genes());
    for gene_id in 0..matrix.n_genes() {
        let key = matrix.gene_symbol(gene_id).to_ascii_uppercase();
        symbol_to_id.entry(key).or_insert(gene_id as u32);
    }

    let mut genesets = Vec::with_capacity(entries.len());
    for (id, (axis, mut symbols)) in entries {
        symbols.sort();
        symbols.dedup();
        let mut gene_ids = Vec::with_capacity(symbols.len());
        let mut missing = Vec::new();
        for symbol in symbols {
            let key = symbol.to_ascii_uppercase();
            if let Some(&gid) = symbol_to_id.get(&key) {
                gene_ids.push(gid);
            } else {
                missing.push(symbol);
            }
        }
        // Sorted by gene index → enables sparse merge in panel_log1p_sum.
        gene_ids.sort_unstable();
        gene_ids.dedup();
        debug!(
            geneset_id = id.as_str(),
            resolved = gene_ids.len(),
            "geneset resolved"
        );
        if !missing.is_empty() {
            let preview: Vec<_> = missing.iter().take(5).cloned().collect();
            warn!(geneset_id = id.as_str(), missing = ?preview, missing_count = missing.len(), "missing genes in geneset");
        }
        genesets.push(Geneset {
            id,
            axis,
            gene_ids,
            missing,
        });
    }

    info!(genesets = genesets.len(), "geneset catalog loaded");

    Ok(GenesetCatalog { genesets })
}

fn parse_catalog_str(
    content: &str,
    source: &str,
    entries: &mut BTreeMap<String, (String, Vec<String>)>,
) -> Result<(), InputError> {
    for (i, line) in content.lines().enumerate() {
        parse_catalog_line(i + 1, line, source, entries)?;
    }
    Ok(())
}

fn parse_catalog_lines<I>(
    lines: I,
    source: String,
    entries: &mut BTreeMap<String, (String, Vec<String>)>,
) -> Result<(), InputError>
where
    I: Iterator<Item = Result<String, std::io::Error>>,
{
    for (i, line) in lines.enumerate() {
        let line = line.map_err(|e| InputError::io(source.clone(), e))?;
        parse_catalog_line(i + 1, &line, &source, entries)?;
    }
    Ok(())
}

fn parse_catalog_line(
    line_no: usize,
    line: &str,
    source: &str,
    entries: &mut BTreeMap<String, (String, Vec<String>)>,
) -> Result<(), InputError> {
    let mut trimmed = line.trim();
    if line_no == 1 {
        trimmed = trimmed.trim_start_matches('\u{feff}');
    }
    if trimmed.is_empty() || trimmed.starts_with('#') {
        return Ok(());
    }
    let cols: Vec<_> = trimmed.split('\t').collect();
    if cols.len() != 3 {
        return Err(InputError::InvalidGenesetCatalog(format!(
            "{source}:{line_no}"
        )));
    }
    if cols[0] == "geneset_id" && cols[1] == "axis" && cols[2] == "gene_symbol" {
        return Ok(());
    }
    let id = cols[0].trim();
    let axis = cols[1].trim();
    let symbol = cols[2].trim();
    if id.is_empty() || axis.is_empty() || symbol.is_empty() {
        return Err(InputError::InvalidGenesetCatalog(format!(
            "{source}:{line_no}"
        )));
    }
    let entry = entries
        .entry(id.to_string())
        .or_insert_with(|| (axis.to_string(), Vec::new()));
    if entry.0 != axis {
        return Err(InputError::InvalidGenesetCatalog(format!(
            "{source}:{line_no}"
        )));
    }
    entry.1.push(symbol.to_string());
    Ok(())
}
