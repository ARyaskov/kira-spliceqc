use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use tracing::{debug, info, warn};

use crate::genesets::aliases::{resolve_entry, symbol_index};

use crate::expression::ExpressionMatrix;
use crate::genesets::{Geneset, GenesetCatalog};
use crate::input::error::InputError;

/// The catalog compiled into the binary, used when no file is found at the
/// catalog path. Exposed so provenance can hash whichever source was used.
pub const EMBEDDED_SPLICE_GENESETS: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/resources/genesets/splicing_genesets.tsv"
));

pub fn load_catalog(
    path: &Path,
    matrix: &dyn ExpressionMatrix,
) -> Result<GenesetCatalog, InputError> {
    // geneset id -> (axis, [(symbol, optional ensembl id)])
    let mut entries: CatalogEntries = BTreeMap::new();
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

    // Case-insensitive index (mouse "Snrpb" resolves) with legacy HGNC
    // aliases (older annotations: "SFRS1" for SRSF1). First occurrence wins.
    let symbol_to_id = symbol_index(matrix);

    let mut genesets = Vec::with_capacity(entries.len());
    let mut aliased_total = 0usize;
    for (id, (axis, mut symbols)) in entries {
        symbols.sort();
        symbols.dedup();
        let mut gene_ids = Vec::with_capacity(symbols.len());
        let mut missing = Vec::new();
        for (symbol, ensembl) in symbols {
            match resolve_entry(&symbol_to_id, &symbol, ensembl.as_deref()) {
                Some((gid, via_alias)) => {
                    if via_alias {
                        aliased_total += 1;
                        debug!(
                            geneset_id = id.as_str(),
                            symbol = symbol.as_str(),
                            "resolved through a legacy alias"
                        );
                    }
                    gene_ids.push(gid);
                }
                None => missing.push(symbol),
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

    if aliased_total > 0 {
        info!(
            aliased = aliased_total,
            "panel genes resolved through legacy symbol aliases"
        );
    }
    info!(genesets = genesets.len(), "geneset catalog loaded");

    Ok(GenesetCatalog { genesets })
}

type CatalogEntries = BTreeMap<String, (String, Vec<(String, Option<String>)>)>;

fn parse_catalog_str(
    content: &str,
    source: &str,
    entries: &mut CatalogEntries,
) -> Result<(), InputError> {
    for (i, line) in content.lines().enumerate() {
        parse_catalog_line(i + 1, line, source, entries)?;
    }
    Ok(())
}

fn parse_catalog_lines<I>(
    lines: I,
    source: String,
    entries: &mut CatalogEntries,
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

/// Catalog line: `geneset_id<TAB>axis<TAB>gene_symbol[<TAB>ensembl_id]`.
fn parse_catalog_line(
    line_no: usize,
    line: &str,
    source: &str,
    entries: &mut CatalogEntries,
) -> Result<(), InputError> {
    let mut trimmed = line.trim();
    if line_no == 1 {
        trimmed = trimmed.trim_start_matches('\u{feff}');
    }
    if trimmed.is_empty() || trimmed.starts_with('#') {
        return Ok(());
    }
    let cols: Vec<_> = trimmed.split('\t').collect();
    if cols.len() != 3 && cols.len() != 4 {
        return Err(InputError::InvalidGenesetCatalog(format!(
            "{source}:{line_no}"
        )));
    }
    if cols[0] == "geneset_id" && cols[1] == "axis" && cols[2] == "gene_symbol" {
        return Ok(()); // header (3 or 4 columns)
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
    let ensembl = cols
        .get(3)
        .map(|c| c.trim())
        .filter(|c| !c.is_empty())
        .map(str::to_string);
    entry.1.push((symbol.to_string(), ensembl));
    Ok(())
}
