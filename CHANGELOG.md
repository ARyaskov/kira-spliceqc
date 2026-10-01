# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

Phases 1 and 3 of the scientific roadmap (Tier A and Tier B, direct
splicing measurements).

### Added

- Input level L2: STARsolo `SJ/` junction count matrices (auto-detected as
  the sibling of a `Gene/` directory or given with `--junctions`), with
  annotation derived from the aligner's `annotated` flag.
- Tier B stage 19: `unannotated_junction_fraction`, `cryptic_3ss_fraction`
  (unannotated acceptors 10-50 nt upstream of an annotated acceptor of the
  same donor, the SF3B1-mutant phenotype), `exon_skip_fraction` (junctions
  skipping an annotated exon vs their inclusion partners) and
  `splice_site_shift` (SpliZ-like), each with a stratum deviation and a
  BH-controlled `_high` flag; `summary.json.junctions` and
  `cells.json.junctions` blocks; `input_levels` reports `L2`.

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
- Control-gene depth correction of every expression panel score (50
  nearest-mean control genes per panel gene, Tirosh et al. 2016) in stages
  2 and 15, and robust standardization within reference stratum and
  library-size bin (`standardize_activity`) for stages 3, 4, 5, 8, 9, 10,
  13 and 15. Null-model library-size correlations of the expression
  signatures drop from 0.4-0.7 to below 0.1; a zero MAD in a bin now gives
  undefined z-scores instead of silent zeros.
- Cell QC flags: `--min-counts` (500) / `--min-genes` (200) mark
  `LOW_DEPTH`, a boolean-like doublet metadata column marks `DOUBLET`;
  both are excluded from reference norms, get undefined deviations and
  appear in `cells.tsv`, `cells.json`, the contract `flags` and
  `summary.json.qc`.
- Provenance block in `summary.json` and `cells.json`: tool version,
  command line, catalog and reference-file CRC-64, input levels, reference
  mode, every model constant and per-metric undefined-cell counts.
- External reference: `kira-spliceqc reference build --input CONTROL --out
  ref.json` stores per-stratum Tier A norms (logit median and
  overdispersion of the unspliced fraction, median and overdispersion of
  the intron retention index, pooled per-gene unspliced ratios);
  `run --reference ref.json` assigns cells to the reference strata and
  computes Tier A deviations and flags against them, so a whole stratum
  shifted relative to the control is detected.
- Legacy HGNC alias table for panel symbol resolution (SR proteins,
  hnRNPs, snRNP/SF3 subunits, U2AF, NMD factors, cell-cycle genes) shared
  by the geneset loader, the stage-15 panels and the cell-cycle scores;
  `summary.json.input.species` and the contract `species` column are
  inferred from symbol casing instead of always `unknown`.
- Stage 18 cell-cycle annotation: Tirosh et al. 2016 S / G2M scores
  (control-gene corrected), Seurat-rule `cell_cycle_phase`, `cycling`
  flag and the `CYCLING` contract flag; `summary.json.cell_cycle` and
  `cells.json.cell_cycle` blocks. Cell-cycle genes are excluded from the
  control pool and `TOP2A` (a G2/M marker) was removed from the
  conflict-risk panel.
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
