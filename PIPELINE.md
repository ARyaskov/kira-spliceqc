# kira-spliceqc Pipeline

This document describes the actual runtime flow and artifacts produced by the current code.

## CLI execution modes

Relevant flags:
- `--run-mode standalone|pipeline` (default: `standalone`)
- `--mode cell|sample` (currently implemented: only `cell`; `sample` returns an error)
- `--extended` enables stages 8-13 (stage 14 is currently not wired in runtime context)
- `--json`, `--tsv` control stage-7 outputs (`both` by default if neither flag is passed)
- `--experimental-signatures` includes composite signatures in stage-7 outputs (implied by `--run-mode pipeline`; `pipeline_step.json` carries `signature_status: "experimental"` and `contract_version`)

## Input resolution (Stage 0)

Order of resolution:

1. If `--run-mode pipeline` and `--cache` provided:
- use the cache file directly (validated via shared-cache dimension checks).

2. If `--run-mode pipeline` and `--input` is a directory:
- detect dataset prefix (if any) and resolve expected shared cache filename:
  - default: `kira-organelle.bin`
  - prefixed dataset: `<PREFIX>.kira-organelle.bin`
- if cache exists and is valid: use shared cache input mode.
- if cache file is missing: log warning and fall back to regular input discovery.
- if cache exists but invalid: hard error (no fallback).

3. Fallback regular detection:
- directory => 10x/MatrixMarket layout via `kira_scio::discover`
- file with `.h5ad` extension => H5AD input

Shared cache spec: [kira-shared-sc-cache/CACHE_FILE.md](https://github.com/ARyaskov/kira-shared-sc-cache/blob/main/CACHE_FILE.md)

## Spliced/unspliced layers (input level L1)

Stage 0 records a layer source in `InputDescriptor.layers`:

1. `--layers PATH` if given (directory with `spliced.mtx`/`unspliced.mtx`, or an `.h5ad` with `layers/`).
2. `spliced.mtx[.gz]` + `unspliced.mtx[.gz]` inside the input directory.
3. STARsolo sibling: `<root>/Gene/<subset>` -> `<root>/Velocyto/<subset>` (also `GeneFull`).
4. `.h5ad`: `layers/spliced` + `layers/unspliced` (+ optional `layers/ambiguous`).

Stage 1 reads the layers in the raw index space, reindexes them with the same
gene/cell sort as the main matrix and keeps them in memory as CSC
(`expression::layers::SplicedUnspliced`). A layer directory with its own
`barcodes.tsv` is matched by barcode (cells absent from it get zero layer counts
and are counted in `cells_without_layers`); otherwise the layers must have the
main matrix's dimensions and order. Gene count mismatches are hard errors.

## Stage order

Runtime order in `run_pipeline`:

- Stage 0: input detection/validation
- Stage 1: expression matrix materialization/opening
- Stage 2: geneset activity aggregation
- Stage 3: isoform entropy/dispersion
- Stage 4: missplicing metrics
- Stage 5: spliceosome imbalance metrics
- Stage 6: SIS (splice integrity score)
- Stage 15: splicing instability proxies (SOS/RLR/SII, expression-only mode)
- Stage 16: Tier A unspliced fraction (only when spliced/unspliced layers were loaded)
- Stages 8-13: only when `--extended`
  - 8 coupling stress
  - 9 exon/intron bias
  - 10 assembly phase imbalance
  - 11 splicing noise
  - 12 cryptic risk
  - 13 spliceosome collapse
- Stage 14: currently not populated in pipeline context (logged as skipped when `--extended`)
- Stage 7: final output serialization
- Pipeline contract generation (`summary.json`, `pipeline_step.json`, `panels_report.tsv`, contract `spliceqc.tsv`) only in `--run-mode pipeline`

## Output directories and artifacts

### Standalone mode

`--out <DIR>` is used directly.

Outputs from stage 7:
- `cells.json` (if `--json` or no explicit format flags)
- `cells.tsv` (if `--tsv` or no explicit format flags)
- `cells.tsv` columns (see the naming convention in METRICS.md):
  - experimental composites: `sis`, `class`, `p_*`
  - expression signatures: `regulator_entropy_expr`, `regulator_dispersion_expr`, `missplicing_burden_expr`, `spliceosome_imbalance_expr`, `coupling_stress_expr`, `exon_definition_bias_expr`, `ea_phase_imbalance_expr`, `b_phase_imbalance_expr`, `catalytic_phase_imbalance_expr`
  - panel cores: `spliceosome_core_expr`, `splicing_rbp_expr`, `rloop_resolution_expr`, `conflict_risk_expr`, `nmd_factor_expr`
  - Tier A (empty without layers): `spliced_umis`, `unspliced_umis`, `ambiguous_umis`, `unspliced_fraction`, `unspliced_fraction_ci_low`, `unspliced_fraction_ci_high`
  - experimental scores: `SOS`, `RLR`, `SII`
  - experimental flags: `splice_overload_high`, `rloop_risk_high`, `splicing_instability_high`, `genome_instability_splicing_flag`

### Pipeline mode

Effective output directory:
- `<OUT>/kira-spliceqc`

Outputs:
- stage-7 outputs (`cells.json`, `cells.tsv`) with same flag rules as standalone
- pipeline contract outputs:
  - `spliceqc.tsv` (contract-formatted table for pipeline integration)
  - `panels_report.tsv`
- `summary.json`
  - includes additive `splicing_instability` block: `panel_version`, `thresholds`, `global_stats`, `cluster_stats`, `missingness`
- `pipeline_step.json`

Note: the contract table `spliceqc.tsv` and the per-cell table `cells.tsv` have distinct names, so neither overwrites the other.

## Stage 1 internal cache: `expr.bin`

`expr.bin` is an internal stage-1 binary cache used when input comes from 10x/H5AD.
It is written to `<effective out dir>/.kira-spliceqc-cache/expr.bin` and the
directory is removed when the run finishes; it is not an output artifact.
When stage-0 selected shared-cache input (`kira-organelle.bin`), stage-1 opens that shared cache directly and does not create `expr.bin`.

Format (little-endian, mmap-friendly):

Header (52 bytes):
- `magic`: `[u8; 8] = b"KIRAEXP1"`
- `version`: `u32 = 1`
- `n_genes`: `u32`
- `n_cells`: `u32`
- `counts_offset`: `u64`
- `libsize_offset`: `u64`
- `gene_index_offset`: `u64`
- `cell_index_offset`: `u64`

Counts section (gene-major sparse rows):
- per gene: `row_offset: u64`, `nnz: u32`, then `nnz` pairs `(cell_id: u32, count: u32)`
- pairs sorted by `cell_id`

Indexes:
- gene symbols and cell names stored as UTF-8 null-terminated strings

Determinism:
- genes sorted lexicographically (stable tie-break by original index)
- cells sorted lexicographically (stable tie-break by original index)
