# kira-spliceqc Pipeline

This document describes the actual runtime flow and artifacts produced by the current code.

## CLI execution modes

Relevant flags:
- `--run-mode standalone|pipeline` (default: `standalone`)
- `--mode cell|sample` (currently implemented: only `cell`; `sample` returns an error)
- `--extended` enables stages 8-13 (stage 14 is currently not wired in runtime context)
- `--json`, `--tsv` control stage-7 outputs (`both` by default if neither flag is passed)

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

## Stage order

Runtime order in `run_pipeline`:

- Stage 0: input detection/validation
- Stage 1: expression matrix materialization/opening
- Stage 2: geneset activity aggregation
- Stage 3: isoform entropy/dispersion
- Stage 4: missplicing metrics
- Stage 5: spliceosome imbalance metrics
- Stage 6: SIS (splice integrity score)
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
- `spliceqc.json` (if `--json` or no explicit format flags)
- `spliceqc.tsv` (if `--tsv` or no explicit format flags)

### Pipeline mode

Effective output directory:
- `<OUT>/kira-spliceqc`

Outputs:
- stage-7 outputs (`spliceqc.json`, `spliceqc.tsv`) with same flag rules as standalone
- pipeline contract outputs:
  - `spliceqc.tsv` (contract-formatted table for pipeline integration)
  - `panels_report.tsv`
  - `summary.json`
  - `pipeline_step.json`

Note: in pipeline mode, contract `spliceqc.tsv` is always written and overwrites any stage-7 TSV with the same name.

## Stage 1 internal cache: `expr.bin`

`expr.bin` is an internal stage-1 binary cache used when input comes from 10x/H5AD.  
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
