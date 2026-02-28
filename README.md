# kira-spliceqc

Deterministic splicing quality control for single-cell RNA-seq.

## Build requirements

- Rust >= 1.95

## Install

Install from crates.io:

```bash
cargo install kira-spliceqc
```

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

## Modes

- `--mode cell` (default): per-cell QC run.
- `--mode sample`: currently not implemented (returns an error).
- `--run-mode standalone` (default): writes stage outputs to `--out`.
- `--run-mode pipeline`: writes into `<OUT>/kira-spliceqc` and generates pipeline contract artifacts.
- `--extended`: enables stages 8-13 (`coupling`, `exon/intron`, `assembly`, `noise`, `cryptic risk`, `collapse`).

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

- `spliceqc.json` (when `--json` is set, or by default when no format flags are passed)
- `spliceqc.tsv` (when `--tsv` is set, or by default when no format flags are passed)

Pipeline mode (`--run-mode pipeline`), output directory: `<OUT>/kira-spliceqc`:

- `spliceqc.tsv` (pipeline contract table)
- `summary.json` (aggregate distributions/regimes/QC fractions)
- `panels_report.tsv` (panel coverage/sum quantiles)
- `pipeline_step.json` (pipeline ingestion manifest)
- `spliceqc.json` (stage-7 JSON output, depending on `--json`/`--tsv` flags)

## Shared cache specification

- Cache format specification: [kira-shared-sc-cache/CACHE_FILE.md](https://github.com/ARyaskov/kira-shared-sc-cache/blob/main/CACHE_FILE.md)
- The shared cache is validated before use (dimensions and format checks from the shared cache reader).

## SIMD note

- SIMD backend is selected at compile time.
- Selected backend is logged at startup.
- Scalar fallback is always available.
