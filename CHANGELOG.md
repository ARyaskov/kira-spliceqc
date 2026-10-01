# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/).

Compatibility policy (full table in `docs/src/compatibility.md`): output
column names and meanings, flag semantics, the pipeline contract, the
geneset catalog and the stage-15 panels change only in a major release;
minor releases may append columns, keys, genesets and flags; patch
releases change none of them. Before 1.0 a minor release may still break
compatibility, and every such change is listed under "Changed (breaking)"
with the old and the new name. Outputs are byte-identical for the same
inputs, flags and version.

## [Unreleased]

## [0.5.0] - 2026-10-01

Phases 1-4 of the scientific roadmap: Tier A and Tier B direct splicing
measurements, the validation package, and release / community
integration. Tier 1 (simulation) validation passes with AUROC 1.0 and
flag FPR below 0.5 % for every spiked effect; tiers 2-3 on public data are
specified in `benchmarks/` and not yet executed.

### Changed (breaking)

- Expression signatures (`*_expr`) are control-gene depth-corrected and
  standardized within reference stratum and library-size bin; every value
  differs from 0.3.0 (same names, same ranges). A zero MAD in a bin gives
  an undefined z-score instead of a silent zero.
- Experimental composite flags (`splice_overload_high`, `rloop_risk_high`,
  `splicing_instability_high`) follow the shared outlier rule (stratum
  deviation >= 3, BH-adjusted p < 0.05) instead of fixed cut-offs;
  `summary.json.splicing_instability.thresholds` describes the rule.
- Stage-15 panel version `SPLICEQC_INSTABILITY_PANEL_V2`: `TOP2A` (a G2/M
  marker) left the conflict-risk panel; cell-cycle genes are excluded from
  every control pool.
- Pipeline contract version `0.5`: the `flags` vocabulary gains
  `LOW_DEPTH`, `DOUBLET` and `CYCLING`; `species` is inferred instead of
  always `unknown`. Column names and order of `spliceqc.tsv` are unchanged.
- Expression cache `expr.bin` format version 2 (gene ids); caches written
  by 0.3.0 are rebuilt automatically.

### Added (release and integration)

- Release workflow: static binaries for Linux (x86_64, aarch64), macOS
  (arm64, x86_64) and Windows on `v*` tags, with a simulate -> run ->
  validate smoke test, SHA-256 sums, release notes from this changelog and
  `cargo publish`.
- `packaging/`: Dockerfile (BioContainers-style), bioconda recipe,
  nf-core-style Nextflow module with an example workflow, and a channel
  overview.
- `python/`: pip-installable wrapper that runs the binary and returns an
  `AnnData` (`cells.tsv` in `obs`, `summary.json` in `uns`).
- `docs/`: mdBook site (installation, inputs, reference model, outputs,
  metric cards with formulas / ranges / confounders / literature, the full
  specification, interpretation rules, "what the tool does not do", three
  tutorials, validation, compatibility policy, release checklist) built to
  GitHub Pages by the docs workflow.
- `resources/genesets/README.md`: catalog version, format and the source
  literature of every geneset; GitHub issue templates for bug reports, new
  panels and new validation datasets.
- `CITATION.cff`, `.zenodo.json` and a JOSS application-note draft
  (`paper/`).
- `kira_spliceqc_mqc.json` in pipeline mode: a MultiQC custom-content table
  (cells, LOW_DEPTH / DOUBLET / cycling fractions, Tier A and Tier B medians
  and flag fractions, reference mode) that `multiqc` renders without a
  dedicated module; listed in `pipeline_step.json.artifacts.multiqc`.

### Added (validation)

- `kira-spliceqc simulate`: Poisson dataset with two cell types, layers,
  junction matrix and spiked effects (cryptic 3' splice sites, intron
  retention, damaged cells, exon skipping) with a `truth.tsv`.
- `kira-spliceqc validate`: AUROC / AUPRC per metric and precision /
  recall / F1 / FPR per flag against a truth table, overall and per cell
  type, as JSON and Markdown; `--pair truth:metric[:flag[:sign]]`.
- `benchmarks/`: tier-1 runner, STARsolo alignment template, dataset
  manifest and acceptance criteria for tiers 2-3; `tests/validation_tier1.rs`
  requires AUROC >= 0.95, flag FPR <= 1 % and recall >= 0.9 for every
  spiked effect.

### Changed

- `exon_skip_fraction` halves the inclusion UMIs (`skip / (skip +
  inclusion / 2)`): an included exon is supported by two junctions and a
  skipped one by a single junction, so the previous `skip / (skip +
  inclusion)` understated skipping by up to a factor of two relative to the
  rMATS junction-count PSI complement it cites.
- `LOW_DEPTH` and `DOUBLET` cells get undefined Tier A and Tier B
  deviations (`unspliced_fraction_dev`, `intron_retention_index_dev`,
  `cryptic_3ss_fraction_dev`, `exon_skip_fraction_dev`,
  `splice_site_shift_dev`) as documented; they were excluded from the
  norms but still received a deviation.
- `reference build` leaves out norms of a stratum without defined cells
  (an empty `global` stratum) instead of writing NaN, which made the file
  unreadable.
- `splice_site_shift_dev` is standardized within stratum and
  junction-depth bin (the raw score rises with junction depth); false
  `splice_site_shift_high` calls on the simulation drop from 7 % to 1 %.
- Experimental composite flags (`splice_overload_high`, `rloop_risk_high`,
  `splicing_instability_high`) use the shared outlier rule (signed
  composite standardized within the stratum, deviation >= 3, BH-adjusted
  p < 0.05) instead of fixed cut-offs that fired on 4 % of null cells;
  the JSON `thresholds` block now states the rule. On the null model the
  three flags fire on 0.0-0.1 % of cells.

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
  `cells.json.junctions` blocks; `input_levels` reports `L2`. The
  null-model test carries a structure-free junction matrix: all three
  Tier B flags fire on <= 0.4 % of null cells.

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
  shifted relative to the control is detected. The file also stores
  per-stratum depth-binned median / MAD norms of every catalog geneset's
  raw activity, of the regulator entropy and of the stage-15 cores, and a
  run with `--reference` standardizes the expression signatures against
  them (`expression_signatures` in `summary.json.reference.external_metrics`,
  `provenance.reference.external_expression_norms`); a file without them
  leaves the signatures dataset-relative.
- Gene ids (Ensembl) are carried through the expression cache (format
  version 2) and indexed next to the symbols; catalogs may carry a fourth
  `ensembl_id` column and are selected with `--catalog`; the species is
  inferred from `ENSG` / `ENSMUSG` / `ENSRNOG` prefixes when ids exist.
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
