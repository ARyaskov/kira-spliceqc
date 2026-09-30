# kira-spliceqc

Deterministic splicing quality control for single-cell RNA-seq.

## Build requirements

- Rust >= 1.95

## Install

Install from crates.io:

```bash
cargo install kira-spliceqc
```

Default geneset catalog is embedded at compile time and used automatically when
`resources/genesets/splicing_genesets.tsv` is unavailable at runtime.

## Usage examples

Standalone run (cell mode, default):

```bash
kira-spliceqc run \
  --input ./data/pbmc3k \
  --out ./out/pbmc3k \
  --mode cell
```

Standalone run with extended stages (8-13):

```bash
kira-spliceqc run \
  --input ./data/pbmc3k \
  --out ./out/pbmc3k \
  --extended
```

Pipeline run (shared cache lookup + pipeline artifacts):

```bash
kira-spliceqc run \
  --input ./data/inf \
  --out ./out/inf \
  --run-mode pipeline \
  --extended
```

Pipeline run with explicit cache path override:

```bash
kira-spliceqc run \
  --input ./data/inf \
  --cache ./data/inf/kira-organelle.bin \
  --out ./out/inf \
  --run-mode pipeline
```

## Input levels

- **L0 — gene counts**: 10x MatrixMarket directory, `.h5ad`, or the shared
  `kira-organelle.bin` cache. Enables the expression signatures (`*_expr`).
- **L1 — spliced / unspliced (/ ambiguous) counts**: auto-detected as
  `spliced.mtx[.gz]` + `unspliced.mtx[.gz]` next to `matrix.mtx` (kb-python),
  in the STARsolo sibling `Solo.out/Velocyto/<subset>` of `Solo.out/Gene/<subset>`,
  or as `layers/spliced` + `layers/unspliced` inside an `.h5ad`. Pass `--layers PATH`
  to point at another directory or file. When the layer directory carries its own
  `barcodes.tsv`, cells are matched by barcode; otherwise the layers must have the
  main matrix's shape and order.

## Modes

- `--mode cell` (default): per-cell QC run.
- `--mode sample`: currently not implemented (returns an error).
- `--run-mode standalone` (default): writes stage outputs to `--out`.
- `--run-mode pipeline`: writes into `<OUT>/kira-spliceqc` and generates pipeline contract artifacts.
- `--layers PATH`: explicit spliced/unspliced layer source (see "Input levels").
- `--metadata PATH`: cell metadata table (`barcode` + columns); auto-detected as `metadata.tsv[.gz]` next to a 10x directory, `.h5ad` inputs use `obs`.
- `--stratify-by COLUMN`: metadata column defining reference strata (default: `cell_type`-like, then `cluster`-like columns, else one global stratum). Strata under 50 cells fold into `global`.
- `--extended`: enables stages 8-13 (`coupling`, `exon/intron`, `assembly`, `noise`, `cryptic risk`, `collapse`).
- `--experimental-signatures`: writes the experimental composite signatures (`sis`/`class`, `SOS`/`RLR`/`SII` and their flags, cryptic risk, collapse) to the per-cell outputs. Off by default; implied by `--run-mode pipeline` because the pipeline contract is built on them.

## Pipeline cache lookup

In pipeline mode, `kira-spliceqc` resolves shared cache according to [kira-shared-sc-cache/CACHE_FILE.md](https://github.com/ARyaskov/kira-shared-sc-cache/blob/main/CACHE_FILE.md):

- no prefix: `kira-organelle.bin`
- prefixed dataset: `<PREFIX>.kira-organelle.bin`

Behavior:

- cache exists and valid: use shared cache input.
- cache missing: warn and fall back to regular input detection (10x/H5AD).
- cache exists but invalid: hard error (no fallback).

`--cache` overrides lookup and uses the provided cache file directly.

## Output artifacts

Standalone mode (`--run-mode standalone`):

- `cells.json` (when `--json` is set, or by default when no format flags are passed)
- `cells.tsv` (when `--tsv` is set, or by default when no format flags are passed)

Pipeline mode (`--run-mode pipeline`), output directory: `<OUT>/kira-spliceqc`:

- `spliceqc.tsv` (pipeline contract table)
- `summary.json` (aggregate distributions/regimes/QC fractions)
- `panels_report.tsv` (panel coverage/sum quantiles)
- `pipeline_step.json` (pipeline ingestion manifest)
- `cells.json` / `cells.tsv` (per-cell stage-7 outputs, depending on `--json`/`--tsv` flags)

Before v0.3 the per-cell stage-7 table was also named `spliceqc.tsv` and was
overwritten by the pipeline contract table in pipeline mode.

## Tier A: unspliced fraction

With spliced/unspliced layers (input level L1) every cell gets
`unspliced_fraction = U / (S + U)` with a Wilson 95 % interval, plus the raw
`spliced_umis`, `unspliced_umis` and `ambiguous_umis`. Fractions are undefined
for cells with fewer than 100 layer UMIs. `unspliced_fraction_dev` is the
logit-scale deviation from the cell's reference stratum and
`nuclear_fraction_flag` marks damaged-cell candidates (fraction far below the
stratum, FDR 5 %). `intron_retention_index` is the per-cell median log2 ratio of
per-gene unspliced ratios to the stratum's pooled ratios (beta-shrunk, 10
pseudo-counts), with `ir_gene_dispersion`, `ir_genes_used`,
`intron_retention_index_dev` and an `intron_retention_high` flag. See METRICS.md
for the reference model and interpretation caveats (protocol and cell-type
dependence).

## Splicing instability proxies

`kira-spliceqc` now emits additive, single-sample-compatible transcriptional proxies for genome/nuclear instability interpretation:

- `SOS` (Spliceosome Overload Score)
- `RLR` (R-loop Risk Proxy)
- `SII` (Splicing Instability Index)

These are deterministic expression-only metrics (no timepoints, no ML). Per-cell values and flags are appended to stage-7 TSV/JSON outputs, and pipeline `summary.json` includes a `splicing_instability` block with thresholds, robust z-score references, quantiles, and missingness.

## Gene symbols and species

Panels are written with current HGNC symbols. Matching is case-insensitive
(mouse `Srsf1` resolves) and falls back to a table of legacy aliases
(`SFRS1` -> `SRSF1`, `ASCC3L1` -> `SNRNP200`, `U2AF65` -> `U2AF2`, ...), so
datasets on older annotations do not silently lose panels. The species in
`summary.json` is inferred from symbol casing (`human` / `mouse` / `unknown`).
Ensembl-ID matching is not yet available.

## Depth correction and reference strata

Panel scores subtract a control-gene background (50 genes of matching mean
expression per panel gene, Tirosh et al. 2016) and are standardized within the
cell's reference stratum and library-size bin. On a Poisson null model this
removes the library-size correlation of the expression signatures (|rho| < 0.1)
and keeps production flags at or below 1 % of cells. `tests/null_model.rs`
enforces both.

## Cell-cycle annotation

Every cell gets Tirosh et al. 2016 S and G2/M scores (control-gene corrected),
a Seurat-rule `cell_cycle_phase` and a `cycling` flag (`CYCLING` in the pipeline
contract). Spliceosome and R-loop panels are enriched for genes that rise in
S/G2M, so cycling cells' expression signatures should be read with that in mind.

## Metric naming

Metrics derived purely from panel expression carry the `_expr` suffix
(`spliceosome_imbalance_expr`, `nmd_factor_expr`, ...). They are expression
signatures, not measurements of splicing. Composite indices (`sis`, `SOS`,
`RLR`, `SII`, regimes) are experimental until validated; see METRICS.md for the
full mapping from v0.2 names.

## Shared cache specification

- Cache format specification: [kira-shared-sc-cache/CACHE_FILE.md](https://github.com/ARyaskov/kira-shared-sc-cache/blob/main/CACHE_FILE.md)
- The shared cache is validated before use (dimensions and format checks from the shared cache reader).

## SIMD note

- SIMD backend is selected at compile time.
- Selected backend is logged at startup.
- Scalar fallback is always available.
