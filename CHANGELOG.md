# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

Phase 1 of the scientific roadmap (Tier A, direct splicing measurements).

### Added

- Input level L1: spliced/unspliced/ambiguous count layers are auto-detected
  next to a 10x directory, in the STARsolo `Velocyto/` sibling of `Gene/`,
  or inside an `.h5ad` (`layers/`), and can be pointed at with `--layers`.
  Layers are reindexed with the main matrix (`SplicedUnspliced`).
- Tier A stage 16: per-cell `unspliced_fraction` with a Wilson 95 % interval
  and raw spliced/unspliced/ambiguous UMI totals in `cells.tsv`/`cells.json`;
  `summary.json` gains `input.levels` and an `unspliced` block. Undefined
  below 100 layer UMIs.
- Cell metadata: `metadata.tsv[.gz]` next to a 10x directory (or
  `--metadata`) and AnnData `obs` string/categorical columns are aligned to
  the canonical cell order.
- Reference strata (`reference` module): `--stratify-by COLUMN`, else
  cell-type / cluster aliases, else global; strata under 50 cells fold into
  `global`. Deviations are robust z-scores (continuous metrics) or
  logit-scale deviations with method-of-moments overdispersion
  (proportions); outlier flags require |d| >= 3 and BH-adjusted p < 0.05.
- Tier A: `unspliced_fraction_dev` and `nuclear_fraction_flag`
  (damaged-cell candidate, DropletQC-style) per cell; per-stratum reference
  medians in `summary.json` and `cells.json`.
- Tier A stage 17: `intron_retention_index` (median log2 ratio of per-gene
  unspliced ratios to the stratum's pooled ratio, beta-shrunk with 10
  pseudo-counts), `ir_gene_dispersion`, `ir_genes_used`,
  `intron_retention_index_dev` and `intron_retention_high`; summary block
  `intron_retention`.

## [0.3.0] - 2026-10-01

Phase 0 of the scientific roadmap: correctness fixes, honest naming and a
null-model guard. Expression signatures are still uncorrected for library
size and cell type; see the Phase 1 targets in `tests/null_model.rs`.

### Changed (breaking)

- Per-cell stage-7 outputs are now `cells.tsv` / `cells.json` in both run
  modes. `spliceqc.tsv` is only the pipeline contract table, so the two no
  longer overwrite each other in pipeline mode.
- Metrics derived purely from panel expression carry the `_expr` suffix
  (`regulator_entropy_expr`, `spliceosome_imbalance_expr`,
  `coupling_stress_expr`, `spliceosome_core_expr`, `nmd_factor_expr`, ...).
  `cells.json` `schema_version` is `2.0`; METRICS.md carries the mapping.
- Composite signatures (`sis`/`class`/penalties, `SOS`/`RLR`/`SII` and their
  flags, cryptic risk, collapse) are experimental and are written only with
  `--experimental-signatures`. Pipeline mode implies the flag; `summary.json`
  and `pipeline_step.json` now say so (`experimental_signatures`,
  `signature_status`, `contract_version = 0.3`).
- Pipeline contract metrics use one strictly monotone map each
  (`unit_clamp`, `saturate`, `rational_saturate`, `sigmoid`) instead of the
  piecewise `normalize01`, which was non-monotone (1.0 -> 1.0 but
  1.2 -> 0.77). Values change; ordering does not.
- `cryptic_risk` is the mean of its three saturated components and spans
  `[0, 1]` (previously confined to `[0.18, 0.82]`).

### Fixed

- Missing contract metrics are written as empty fields with a
  `MISSING_METRICS` flag and `regime = Unclassified`; they are no longer
  coerced to `0`, which had read as "worst fidelity, no stress".
- Stage 3 used a private `1e-12` epsilon and no zero-MAD guard for the
  entropy z-score, so a collapsed MAD produced z ~ 1e12 and saturated the
  SIS penalty. It now uses the shared robust z-score.
- Stage 4 accepted two resolved core panels but produced NaN burden for
  every cell unless all three were present. The core mean now averages the
  finite core panels (minimum two).
- The stage-1 expression cache (`expr.bin`) no longer remains in the output
  directory.

### Added

- Warning when a panel's MAD is zero (its z-scores collapse to 0 and it
  carries no signal).
- `tests/null_model.rs`: deterministic Poisson dataset with no splicing
  biology; asserts determinism, absence of missing-value floods, and Phase 0
  baselines for flag fractions and library-size correlation.
- GitHub Actions CI (build, clippy `-D warnings` on lib/bin, tests on Linux
  and macOS, null-model report).
- `CHANGELOG.md`.

## [0.2.0]

- Additional metrics (stages 8-15) and pipeline contract.

## [0.1.0]

- Initial release.
