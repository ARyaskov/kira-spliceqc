use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use tracing::{info, warn};

use crate::expression::ExpressionMatrix;
use crate::genesets::{Geneset, GenesetCatalog};
use crate::input::error::InputError;

pub fn load_catalog(
    path: &Path,
    matrix: &dyn ExpressionMatrix,
) -> Result<GenesetCatalog, InputError> {
    let file = File::open(path)
        .map_err(|_| InputError::GenesetCatalogMissing(path.display().to_string()))?;
    let reader = BufReader::new(file);

    let mut entries: BTreeMap<String, (String, Vec<String>)> = BTreeMap::new();
    let mut line_no = 0usize;
    for line in reader.lines() {
        line_no += 1;
        let line = line.map_err(|e| InputError::io(path, e))?;
        let mut trimmed = line.trim();
        if line_no == 1 {
            trimmed = trimmed.trim_start_matches('\u{feff}');
        }
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let cols: Vec<_> = trimmed.split('\t').collect();
        if cols.len() != 3 {
            return Err(InputError::InvalidGenesetCatalog(format!(
                "{}:{}",
                path.display(),
                line_no
            )));
        }
        if cols[0] == "geneset_id" && cols[1] == "axis" && cols[2] == "gene_symbol" {
            continue;
        }
        let id = cols[0].trim();
        let axis = cols[1].trim();
        let symbol = cols[2].trim();
        if id.is_empty() || axis.is_empty() || symbol.is_empty() {
            return Err(InputError::InvalidGenesetCatalog(format!(
                "{}:{}",
                path.display(),
                line_no
            )));
        }
        let entry = entries
            .entry(id.to_string())
            .or_insert_with(|| (axis.to_string(), Vec::new()));
        if entry.0 != axis {
            return Err(InputError::InvalidGenesetCatalog(format!(
                "{}:{}",
                path.display(),
                line_no
            )));
        }
        entry.1.push(symbol.to_string());
    }

    let mut symbol_to_id = std::collections::HashMap::new();
    for gene_id in 0..matrix.n_genes() {
        symbol_to_id.insert(matrix.gene_symbol(gene_id).to_string(), gene_id as u32);
    }

    let mut genesets = Vec::with_capacity(entries.len());
    for (id, (axis, mut symbols)) in entries {
        symbols.sort();
        let mut gene_ids = Vec::new();
        let mut missing = Vec::new();
        for symbol in symbols {
            if let Some(&gid) = symbol_to_id.get(&symbol) {
                gene_ids.push(gid);
            } else {
                missing.push(symbol);
            }
        }
        info!(
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
