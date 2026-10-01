//! MultiQC custom-content export of the pipeline summary.
//!
//! `kira_spliceqc_mqc.json` is one row per sample in MultiQC's
//! [custom content](https://docs.seqera.io/multiqc/custom_content) table
//! format: `multiqc <dir>` picks it up without a dedicated module. Only
//! sample-level numbers are exported (per-cell tables stay in `cells.tsv`);
//! Tier A / Tier B columns are present only when their input level was
//! available.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use serde_json::{Map, Value, json};

use crate::input::error::InputError;

pub const MULTIQC_FILE: &str = "kira_spliceqc_mqc.json";

/// Column order and headers of the MultiQC table.
const COLUMNS: &[(&str, &str, &str, &str)] = &[
    // key, title, description, format
    ("n_cells", "Cells", "Cells in the input matrix", "{:,.0f}"),
    (
        "low_depth_fraction",
        "Low depth",
        "Fraction of cells flagged LOW_DEPTH (below --min-counts / --min-genes)",
        "{:.1%}",
    ),
    (
        "doublet_fraction",
        "Doublets",
        "Fraction of cells flagged DOUBLET from metadata",
        "{:.1%}",
    ),
    (
        "cycling_fraction",
        "Cycling",
        "Fraction of cells in S or G2M phase",
        "{:.1%}",
    ),
    (
        "unspliced_fraction_median",
        "Unspliced",
        "Median unspliced fraction U / (S + U) (Tier A)",
        "{:.3f}",
    ),
    (
        "nuclear_fraction_flag_fraction",
        "Damaged",
        "Fraction of cells with nuclear_fraction_flag (Tier A)",
        "{:.1%}",
    ),
    (
        "intron_retention_index_median",
        "IR index",
        "Median intron retention index, log2 vs stratum (Tier A)",
        "{:.3f}",
    ),
    (
        "intron_retention_high_fraction",
        "IR high",
        "Fraction of cells with intron_retention_high (Tier A)",
        "{:.1%}",
    ),
    (
        "cryptic_3ss_fraction_median",
        "Cryptic 3'SS",
        "Median cryptic 3' splice-site usage (Tier B)",
        "{:.4f}",
    ),
    (
        "cryptic_3ss_high_fraction",
        "Cryptic high",
        "Fraction of cells with cryptic_3ss_high (Tier B)",
        "{:.1%}",
    ),
    (
        "exon_skip_fraction_median",
        "Exon skip",
        "Median exon skip fraction (Tier B)",
        "{:.4f}",
    ),
    (
        "exon_skip_high_fraction",
        "Skip high",
        "Fraction of cells with exon_skip_high (Tier B)",
        "{:.1%}",
    ),
    (
        "splice_site_shift_high_fraction",
        "Shift high",
        "Fraction of cells with splice_site_shift_high (Tier B)",
        "{:.1%}",
    ),
    (
        "reference_mode",
        "Reference",
        "Reference mode: global, stratified or external",
        "{}",
    ),
];

/// Builds the MultiQC document for one sample from a parsed `summary.json`.
pub fn multiqc_document(sample: &str, summary: &Value) -> Value {
    let mut row = Map::new();
    let mut put = |key: &str, value: Option<Value>| {
        if let Some(v) = value.filter(|v| !v.is_null()) {
            row.insert(key.to_string(), v);
        }
    };
    put("n_cells", summary.pointer("/input/n_cells").cloned());
    put(
        "low_depth_fraction",
        summary.pointer("/qc/low_depth_fraction").cloned(),
    );
    put(
        "doublet_fraction",
        summary.pointer("/qc/doublet_fraction").cloned(),
    );
    put(
        "cycling_fraction",
        summary.pointer("/cell_cycle/cycling_fraction").cloned(),
    );
    put(
        "unspliced_fraction_median",
        summary.pointer("/unspliced/median").cloned(),
    );
    put(
        "nuclear_fraction_flag_fraction",
        summary
            .pointer("/unspliced/nuclear_fraction_flag_fraction")
            .cloned(),
    );
    put(
        "intron_retention_index_median",
        summary.pointer("/intron_retention/median").cloned(),
    );
    put(
        "intron_retention_high_fraction",
        summary
            .pointer("/intron_retention/high_flag_fraction")
            .cloned(),
    );
    put(
        "cryptic_3ss_fraction_median",
        summary
            .pointer("/junctions/cryptic_3ss_fraction_median")
            .cloned(),
    );
    put(
        "cryptic_3ss_high_fraction",
        summary
            .pointer("/junctions/cryptic_3ss_high_fraction")
            .cloned(),
    );
    put(
        "exon_skip_fraction_median",
        summary
            .pointer("/junctions/exon_skip_fraction_median")
            .cloned(),
    );
    put(
        "exon_skip_high_fraction",
        summary
            .pointer("/junctions/exon_skip_high_fraction")
            .cloned(),
    );
    put(
        "splice_site_shift_high_fraction",
        summary
            .pointer("/junctions/splice_site_shift_high_fraction")
            .cloned(),
    );
    put(
        "reference_mode",
        summary.pointer("/reference/mode").cloned(),
    );

    let mut headers = Map::new();
    for (key, title, description, format) in COLUMNS {
        if !row.contains_key(*key) {
            continue;
        }
        let mut header = Map::new();
        header.insert("title".into(), json!(title));
        header.insert("description".into(), json!(description));
        header.insert("format".into(), json!(format));
        if format.ends_with("%}") {
            header.insert("min".into(), json!(0));
            header.insert("max".into(), json!(1));
            header.insert("scale".into(), json!("OrRd"));
        }
        headers.insert((*key).to_string(), Value::Object(header));
    }

    let version = summary
        .pointer("/tool/version")
        .and_then(Value::as_str)
        .unwrap_or(env!("CARGO_PKG_VERSION"));
    json!({
        "id": "kira_spliceqc",
        "section_name": "kira-spliceqc",
        "section_href": "https://github.com/ARyaskov/kira-spliceqc",
        "description": format!(
            "Splicing QC summary (kira-spliceqc {version}). Fractions are shares of cells; \
             deviations and flags are relative to the reference strata recorded in summary.json."
        ),
        "plot_type": "table",
        "pconfig": {
            "id": "kira_spliceqc_table",
            "title": "kira-spliceqc: splicing QC",
            "namespace": "kira-spliceqc",
        },
        "headers": headers,
        "data": { sample: row },
    })
}

/// Writes `kira_spliceqc_mqc.json` next to `summary.json`.
pub fn write_multiqc(out_dir: &Path, sample: &str, summary: &Value) -> Result<(), InputError> {
    let path = out_dir.join(MULTIQC_FILE);
    let file = File::create(&path).map_err(|e| InputError::io(&path, e))?;
    let mut w = BufWriter::with_capacity(1 << 14, file);
    serde_json::to_writer_pretty(&mut w, &multiqc_document(sample, summary))
        .map_err(|e| InputError::OutputSerialization(e.to_string()))?;
    w.write_all(b"\n").map_err(|e| InputError::io(&path, e))?;
    w.flush().map_err(|e| InputError::io(&path, e))?;
    Ok(())
}

/// Sample name for the MultiQC row: the input's file stem (`pbmc3k` for
/// `./data/pbmc3k` or `./data/pbmc3k.h5ad`), `sample` when it has none.
pub fn sample_name(input: &str) -> String {
    let trimmed = input.trim_end_matches(['/', '\\']);
    let path = Path::new(trimmed);
    let name = if path
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("h5ad"))
    {
        path.file_stem()
    } else {
        path.file_name()
    };
    match name
        .and_then(|n| n.to_str())
        .filter(|n| !n.is_empty() && *n != ".")
    {
        Some(n) => n.to_string(),
        None => "sample".to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tier_columns_follow_the_input_levels() {
        let summary = json!({
            "tool": {"version": "0.3.0"},
            "input": {"n_cells": 10},
            "qc": {"low_depth_fraction": 0.1, "doublet_fraction": 0.0},
            "cell_cycle": {"cycling_fraction": 0.2},
            "unspliced": null,
            "intron_retention": null,
            "junctions": null,
            "reference": {"mode": "internal"}
        });
        let doc = multiqc_document("s1", &summary);
        let row = &doc["data"]["s1"];
        assert_eq!(row["n_cells"], 10);
        assert!(row.get("unspliced_fraction_median").is_none());
        assert!(doc["headers"].get("unspliced_fraction_median").is_none());
        assert_eq!(doc["headers"]["low_depth_fraction"]["max"], 1);
        assert_eq!(doc["plot_type"], "table");
    }

    #[test]
    fn sample_names() {
        assert_eq!(sample_name("./data/pbmc3k/"), "pbmc3k");
        assert_eq!(sample_name("data/pbmc3k.h5ad"), "pbmc3k");
        assert_eq!(sample_name("Solo.out/Gene/raw"), "raw");
        assert_eq!(sample_name(""), "sample");
    }
}
