# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/).

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
