# kira-spliceqc Metrics Specification

This document defines the metrics produced by `kira-spliceqc`, including formulas, constants, and classification rules.

Scope:
- standalone outputs: `cells.tsv`, `cells.json`
- pipeline outputs: `kira-spliceqc/spliceqc.tsv`, `kira-spliceqc/summary.json`, `kira-spliceqc/pipeline_step.json`, `kira-spliceqc/panels_report.tsv`

## Canonical Conventions

1. Determinism
- No stochastic steps.
- Fixed formulas and fixed thresholds/constants.

2. Axis semantics
- Core metrics are per cell.
- Timecourse metrics are per timepoint aggregate (medians), then converted to deltas.

3. Missing values
- Unresolved values are represented as `NaN` internally and serialized as empty TSV fields / `null` in JSON.

4. Robust normalization primitive
- Many stages use robust z-score:
  - `median(x)` over finite values only.
  - `MAD(x) = median(|x - median(x)|)` over finite values only.
  - `robust_z(x) = (x - median) / (MAD * 1.4826 + EPS)`.

## Naming Convention (v0.3)

Every per-cell metric that is derived purely from the expression of a gene
panel carries the `_expr` suffix. These are expression signatures: they do not
observe junctions, introns or isoforms and must not be read as splicing
measurements. Composite indices built on top of them (`sis`, `class`, `SOS`,
`RLR`, `SII`, `cryptic_risk`, `collapse_status`, contract `regime`) are
experimental until validated on datasets with known splicing defects.

Column / key mapping from v0.2:

| v0.2 name | v0.3 name (cells.tsv) | cells.json location |
| --- | --- | --- |
| `iso_entropy` | `regulator_entropy_expr` | `regulator_expr.entropy` |
| `iso_dispersion` | `regulator_dispersion_expr` | `regulator_expr.dispersion` |
| `missplicing_burden` | `missplicing_burden_expr` | `missplicing_expr.burden` |
| `imbalance` | `spliceosome_imbalance_expr` | `imbalance_expr.imbalance` |
| `coupling_stress` | `coupling_stress_expr` | `coupling_expr.coupling_stress` |
| `exon_definition_bias` | `exon_definition_bias_expr` | `exon_intron_bias_expr.exon_definition_bias` |
| `ea_imbalance` / `b_imbalance` / `cat_imbalance` | `ea_phase_imbalance_expr` / `b_phase_imbalance_expr` / `catalytic_phase_imbalance_expr` | `assembly_phase_expr.*` |
| `splice_core` | `spliceosome_core_expr` | `splicing_instability.spliceosome_core_expr` |
| `rbp_core` | `splicing_rbp_expr` | `splicing_instability.splicing_rbp_expr` |
| `rloop_resolve_core` | `rloop_resolution_expr` | `splicing_instability.rloop_resolution_expr` |
| `conflict_risk_core` | `conflict_risk_expr` | `splicing_instability.conflict_risk_expr` |
| `nmd_core` | `nmd_factor_expr` | `splicing_instability.nmd_factor_expr` |

The formulas below keep their internal symbol names (`iso_entropy`, `imbalance`,
...) for readability; output columns follow the table above. `cells.json`
`schema_version` is `2.0`.

## Notation

Let:
- `count(g, c)` = raw UMI count for gene `g` in cell `c`.
- `libsize(c)` = sum of counts in cell `c` (clamped with `max(1, libsize)` where needed).
- `cp10k(g, c) = 1e4 * count(g, c) / max(1, libsize(c))`.
- `log_cp10k(g, c) = ln(1 + cp10k(g, c))`.
- `relu(x) = max(x, 0)`.
- `sigmoid(x) = 1 / (1 + exp(-x))`.

## Reference Strata and Deviations

Every `_dev` metric and outlier flag is computed relative to the reference of
the cell's own stratum. `summary.json.reference.mode` records which mode ran:

| Mode | When | Norm |
| --- | --- | --- |
| `stratified` | metadata column found (`--stratify-by`, else the first of `cell_type`, `celltype`, `cell_type_annotation`, `annotation`, ..., then `cluster`, `leiden`, `louvain`, `seurat_clusters`, ...) | per stratum; strata with fewer than `MIN_STRATUM_CELLS = 50` cells (and cells with an empty value) fold into the `global` stratum |
| `global` | no usable column | one stratum, the whole dataset |
| `external` | reserved for `--reference ref.json` (not yet available) | |

Metadata sources: `metadata.tsv[.gz]` next to a 10x directory (header line,
first column = barcode) or `--metadata PATH`; `obs` string and categorical
columns of an AnnData input.

Continuous metric `m`, stratum `s(c)`:
- `robust_z(m, c) = (m_c - median_s(m)) / (1.4826 * MAD_s(m))`; NaN when the stratum MAD is 0.

Proportion `p_c` with `n_c` trials (e.g. unspliced fraction with `S + U` UMIs):
- `p_c` is clamped to `[0.5/n_c, 1 - 0.5/n_c]`
- `d_c = (logit(p_c) - median_s(logit p)) / sqrt(1 / (n_c p_c (1 - p_c)) + tau2_s)`
- `tau2_s = max(0, (1.4826 * MAD_s(logit p))^2 - mean_s(1 / (n p (1 - p))))` (method-of-moments overdispersion)

Outlier flags: `sign * d_c >= 3` **and** Benjamini-Hochberg adjusted two-sided
normal p-value `< 0.05` within the stratum. Expected flag rate on a null model
is therefore below 1 %.

## Tier A: Unspliced Fraction (Stage 16, requires input level L1)

Direct measurement from spliced/unspliced count layers (see README "Input
levels"). Not an expression signature: no `_expr` suffix.

Per cell `c`, with `S_c`, `U_c`, `A_c` the spliced, unspliced and ambiguous
UMI totals:
- `spliced_umis = S_c`, `unspliced_umis = U_c`, `ambiguous_umis = A_c` (0 without an ambiguous layer)
- `unspliced_fraction = U_c / (S_c + U_c)`; ambiguous UMIs are excluded from both terms
- `unspliced_fraction_ci_low/high`: Wilson 95 % interval of `U_c / (S_c + U_c)`
- undefined (empty) when `S_c + U_c < MIN_LAYER_UMIS = 100`
- cells absent from a barcode-matched layer directory have `S = U = 0` and are
  counted in `cells_without_layers`
- `unspliced_fraction_dev = d_c` from "Reference Strata and Deviations" with `n_c = S_c + U_c`
- `nuclear_fraction_flag`: `d_c <= -3` and BH-adjusted p < 0.05 within the stratum, i.e. an
  unspliced fraction far below the cell's own stratum: cytoplasmic debris / damaged-cell
  candidate (Muskovic & Powell 2021, DropletQC). Not raised for high fractions.

Interpretation: the expected level depends on the protocol (single-nucleus
0.5-0.7, whole-cell 3' 0.1-0.3) and on cell type; compare within a protocol
and a reference stratum. 3' libraries also count internal priming on A-rich
introns as unspliced (La Manno et al. 2018; Muskovic & Powell 2021).

`summary.json.unspliced` reports `n_defined_cells`, `median`, `p10`, `p90`
(round((n-1)q) quantiles over defined cells), `nuclear_fraction_flag_fraction`,
the layer source and per-stratum median/MAD.

## Tier A: Intron Retention Index (Stage 17, requires input level L1)

Per gene `g`, cell `c`, reference stratum `s = s(c)`:
- `p_gs` = pooled unspliced ratio `sum_c U_gc / sum_c (S_gc + U_gc)` over the stratum's
  cells with `S_gc + U_gc >= MIN_GENE_UMIS = 5`; undefined unless at least
  `MIN_CELLS_PER_GENE = 10` such cells exist; clamped to `[0.5/N, 1 - 0.5/N]`
- `IR_gc = (U_gc + K p_gs) / (S_gc + U_gc + K p_gs + K (1 - p_gs))` with `K = PRIOR_STRENGTH = 10`
  pseudo-counts (beta shrinkage toward the stratum ratio; fixed and modest so that outlier
  cells keep their signal while low-count genes are stabilised)
- gene set `G_c` = genes with `S_gc + U_gc >= 5` and a defined `p_gs`; the cell is undefined
  when `|G_c| < MIN_GENES = 20`
- `intron_retention_index = median_{g in G_c} log2(IR_gc / p_gs)` (log2 units; 0 = at the
  stratum reference, +1 = twice the reference unspliced ratio)
- `ir_gene_dispersion = MAD_{g in G_c} log2(IR_gc / p_gs)`: small when every gene shifts
  together (global retention), large for gene-specific changes
- `ir_genes_used = |G_c|`
- `intron_retention_index_dev` = robust z-score within the stratum; `intron_retention_high`
  = `dev >= 3` and BH-adjusted p < 0.05

Adapted from the IRFinder intron-retention ratio (Middleton et al. 2017 Genome
Biology) to sparse per-cell counts. 3' libraries confound unspliced signal with
internal priming; the per-gene reference absorbs gene-specific background, so
the index measures departure from the stratum, not absolute retention.

## Geneset Activity (Stage 2)

For each geneset `S` and cell `c`:
- `A(S, c) = mean_{g in S}( ln(1 + 1e4 * count(g, c) / max(1, libsize(c))) )`

This matrix is the source for downstream robust z-score metrics.

## Isoform Dispersion (Stage 3)

Regulator union:
- `R = U(SRSF_SR, HNRNP, U1_CORE, U2_CORE, SF3B_AXIS, MINOR_U12)`

For each cell `c`:
- `sum_R = sum_{g in R} cp10k(g, c)`
- `p_g = cp10k(g, c) / (sum_R + EPS_STAGE3)`
- `iso_entropy = -sum_{g in R}(p_g * ln(p_g + EPS_STAGE3)) / ln(|R|)`
- `iso_dispersion = (1 / (sum_{g in R} p_g^2 + EPS_STAGE3)) / |R|`
- `z_entropy = robust_z(iso_entropy)`

If `sum_R <= 0`, `iso_entropy` and `iso_dispersion` are `NaN`.

## Missplicing Burden (Stage 4)

Required panel ids:
- `U1_CORE`, `U2_CORE`, `SF3B_AXIS`, `SRSF_SR`, `HNRNP`, `MINOR_U12`, `NMD_SURVEILLANCE`

Using per-panel robust z-scores:
- `b_core = relu(-mean(finite values among z_u1, z_u2, z_sf3b))`; NaN if fewer than 2 are finite (matches the stage gate of >= 2 resolved core panels)
- `b_u12 = relu(z_u12 - mean(z_u1, z_u2))`
- `b_nmd = relu(z_nmd)`
- `b_srhn = |z_srsf - z_hnrnp|`
- `missplicing_burden = 0.35*b_core + 0.25*b_u12 + 0.25*b_nmd + 0.15*b_srhn`
- `burden_star (splice_junction_noise primitive) = 1 - exp(-missplicing_burden)`

## Spliceosome Imbalance (Stage 5)

Required panel ids:
- `U1_CORE`, `U2_CORE`, `SF3B_AXIS`, `SRSF_SR`, `HNRNP`, `MINOR_U12`, `NMD_SURVEILLANCE`

Axes:
- `axis_sr_hnrnp = z_srsf - z_hnrnp`
- `axis_u2_u1 = z_u2 - z_u1`
- `axis_u12_major = z_u12 - mean(z_u1, z_u2)`
- `axis_nmd = z_nmd`

Integrated imbalance:
- Clamp each of `z_u1, z_u2, z_sf3b, z_srsf, z_hnrnp, z_u12` into `[-6, 6]`.
- `imbalance = sqrt(mean(clamped_z_i^2))`

## Splice Integrity Score (SIS, Stage 6)

Penalty components:
- `p_missplicing = burden_star`
- `p_imbalance = clamp01((imbalance - 0.8) / 1.2)`
- `p_entropy_z = clamp01((z_entropy - 1.5) / 2.0)`
- `p_entropy_abs = clamp01((iso_entropy - 0.85) / 0.15)`

Score:
- `sis_raw = 1 - (0.35*p_missplicing + 0.25*p_imbalance + 0.25*p_entropy_z + 0.15*p_entropy_abs)`
- `sis = clamp01(sis_raw)`

Class:
- `Intact` if `sis >= 0.80`
- `Stressed` if `0.60 <= sis < 0.80`
- `Impaired` if `0.40 <= sis < 0.60`
- `Broken` if `sis < 0.40` or if SIS inputs are non-finite

## Splicing Instability Proxies (Stage 15, expression-only)

Panel version:
- `SPLICEQC_INSTABILITY_PANEL_V1`

Panels (human symbols, stable order):
- Core spliceosome / snRNP load proxy:
  - `SNRPB,SNRPD1,SNRPD2,SNRPD3,SNRPE,SNRPF,SNRPG,SF3A1,SF3A2,SF3A3,SF3B1,SF3B2,SF3B3,SF3B4,SF3B5,PRPF3,PRPF4,PRPF6,PRPF8,PRPF19,U2AF1,U2AF2`
- Splicing regulation / stress-sensitive RBPs:
  - `HNRNPA1,HNRNPA2B1,HNRNPC,HNRNPK,SRSF1,SRSF2,SRSF3,SRSF6,SRSF7,RBM39,RBM10,RBM17`
- R-loop resolution (protective axis):
  - `SETX,DDX5,DDX21,DHX9,RNASEH1,RNASEH2A,RNASEH2B,RNASEH2C,BRCA1,BRCA2`
- Transcription-replication conflict risk (optional):
  - `TOP1,TOP2A,TOP2B,POLR2A,SUPT5H,SUPT6H`
- NMD surveillance (optional):
  - `UPF1,UPF2,UPF3B,SMG1,SMG5,SMG6,SMG7`

Panel trimmed mean:
- collect `v_i = log_cp10k(g_i, c)` for mapped genes in panel
- if mapped values `< MIN_GENES_PER_PANEL_CELL (3)`, score = `NaN`
- `TM(P,c) = mean(v_sorted[k:n-k])`, `k = floor(0.1*n)`

Robust z-score:
- `Z(x) = (x - median(x)) / (1.4826 * MAD(x) + EPS_ROBUST)`
- if `MAD == 0`, z-score for finite values is `0` and a warning naming the panel is logged (the panel then carries no signal)

Spliceosome Overload Score (SOS):
- `splice_core = TM(SpliceosomePanel, c)`
- `rbp_core = TM(SplicingRBPPanel, c)`
- `SOS(c) = 0.65 * Z(splice_core) + 0.35 * Z(rbp_core)`

R-loop Risk Proxy (RLR):
- `rloop_resolve = TM(RloopResolutionPanel, c)`
- `conflict_risk = TM(ConflictPanel, c)` when conflict panel enabled
- if conflict panel enabled:
  - `RLR(c) = 0.7 * relu(-Z(rloop_resolve)) + 0.3 * relu(Z(conflict_risk))`
- if conflict panel disabled:
  - `RLR(c) = relu(-Z(rloop_resolve))`

Splicing Instability Index (SII):
- with NMD panel enabled:
  - `SII(c) = 0.6 * relu(SOS(c)) + 0.4 * relu(-Z(nmd_core))`
- with NMD panel disabled:
  - `SII(c) = relu(SOS(c))`

Flags:
- `splice_overload_high`: `SOS >= 2.0`
- `rloop_risk_high`: `RLR >= 1.5`
- `splicing_instability_high`: `SII >= 2.0`
- `genome_instability_splicing_flag`: `splice_overload_high && rloop_risk_high`

NaN behavior:
- flags evaluate to `false` for non-finite scores
- missingness counters are emitted in JSON outputs

Optional junction-aware mode (future placeholder):
- `JE` (junction entropy)
- `IRB` (intron retention burden)
- `CSP` (cryptic splicing proxy)

## Coupling Stress (Stage 8)

Required panel ids:
- `TRANSCRIPTION_COUPLING`, `U1_CORE`, `U2_CORE`, `SF3B_AXIS`

Metric:
- `coupling_stress = z_transcription_coupling - mean(z_u1, z_u2, z_sf3b)`

## Exon/Intron Definition Bias (Stage 9)

Required panel ids:
- `SRSF_SR`, `HNRNP`, `U2AF_AXIS`

Metric:
- `exon_definition_bias = (z_srsf + z_u2af) - z_hnrnp`

## Assembly Phase Imbalance (Stage 10)

Required panel ids:
- `SPLICE_EA_PHASE`, `SPLICE_B_PHASE`, `SPLICE_CATALYTIC_PHASE`

Metrics:
- `ea_imbalance = z_ea - mean(z_b, z_cat)`
- `b_imbalance = z_b - mean(z_ea, z_cat)`
- `cat_imbalance = z_cat - mean(z_ea, z_b)`

## Splicing Noise (Stage 11)

Core panel ids:
- `U1_CORE`, `U2_CORE`, `SF3B_AXIS`, `SRSF_SR`, `HNRNP`, `MINOR_U12`

Per-panel noise:
- `noise(panel) = MAD(values_panel) / (|median(values_panel)| + EPS_STAGE11)`

Global index:
- `noise_index = mean(noise(panel))` over finite panel noises.

## Cryptic Splicing Risk (Stage 12)

Inputs:
- `axis_sr_hnrnp`, `z_entropy`, `z_nmd`

Saturation normalization with `SAT = 6`:
- `x_sr_hnrnp = sat01(|axis_sr_hnrnp|)`
- `x_entropy = sat01(z_entropy)`
- `x_nmd = sat01(z_nmd)`
- where `sat01(v) = 0 if v<=0; v/6 if 0<v<6; 1 if v>=6`

Risk:
- `cryptic_risk = (x_sr_hnrnp + x_entropy + x_nmd) / 3` (spans the full `[0, 1]` range;
  the summary flags cells with `cryptic_risk > 0.7`)

## Spliceosome Collapse (Stage 13)

Boolean conditions:
- `core_suppression`: `z_u1 < -1.5 && z_u2 < -1.5 && z_sf3b < -1.5`
- `high_imbalance`: `imbalance > 1.8`
- `low_sis`: `sis < 0.4`

Status:
- `Collapse` if all three conditions are true.
- `NoCollapse` otherwise.
- `Inconclusive` if any of `z_u1/z_u2/z_sf3b` is non-finite.

## Timecourse Coherence (Stage 14)

Given timepoint medians:
- `sis_median[t]`, `entropy_median[t]`, `imbalance_median[t]`

Deltas:
- `delta_sis[t] = sis_median[t+1] - sis_median[t]`
- `delta_entropy[t] = entropy_median[t+1] - entropy_median[t]`
- `delta_imbalance[t] = imbalance_median[t+1] - imbalance_median[t]`

Trajectory class (`n_timepoints >= 3` and all finite):
- `Adaptive` if all `delta_sis > 0` and all `delta_entropy <= 0` and all `delta_imbalance <= 0`
- `Degenerative` if all `delta_sis < 0` and all `delta_entropy >= 0` and all `delta_imbalance >= 0`
- `Oscillatory` otherwise
- `Inconclusive` if preconditions are not met

## Pipeline Contract Derived Metrics

Row metrics in `kira-spliceqc/spliceqc.tsv`. Each is a single strictly monotone
map of one source metric onto `[0, 1]`; non-finite sources stay non-finite and
are written as empty fields (see "Missing values" below).

- `splice_fidelity_index = unit_clamp(sis)`
- `intron_retention_rate = saturate(b_u12)`
- `exon_skipping_rate = saturate(|axis_u2_u1|)`
- `alt_splice_burden = rational_saturate(missplicing_burden)`
- `splice_junction_noise = unit_clamp(burden_star)`
- `stress_splicing_index = sigmoid(coupling_stress)` if stage 8 present, else `saturate(imbalance)`

Transforms:
- `unit_clamp(x)`: identity on `[0,1]`, clamped outside (for metrics that are unit-bounded by construction).
- `saturate(x) = 1 - exp(-max(x, 0))` for non-negative unbounded inputs.
- `rational_saturate(x) = max(x, 0) / (1 + max(x, 0))`.
- `sigmoid(x) = 1 / (1 + exp(-x))` for signed inputs.

Note: `alt_splice_burden` and `splice_junction_noise` are two monotone
transforms of the same `missplicing_burden` and carry no independent
information. Both are kept for contract compatibility and are deprecated.

Missing values:
- If any of the six row metrics is non-finite, the row gets the `MISSING_METRICS` flag,
  `regime = Unclassified`, and `confidence` is empty.
- Missing values are never coerced to `0`.

Confidence:
- `penalty = 0.4*intron_retention_rate + 0.3*alt_splice_burden + 0.3*splice_junction_noise`
- `confidence = clamp01(0.65*splice_fidelity_index + 0.35*(1 - penalty))`

Regime classification:
- `SplicingCollapse` if `splice_fidelity_index < 0.2` and `stress_splicing_index > 0.8`
- `SpliceNoiseDominant` if `splice_junction_noise > 0.75`
- `StressInducedSplicing` if `stress_splicing_index > 0.65`
- `RegulatedAlternativeSplicing` if any of `alt_splice_burden`, `exon_skipping_rate`, `intron_retention_rate` is `> 0.45`
- `HighFidelitySplicing` if `splice_fidelity_index > 0.75` and `splice_junction_noise < 0.35` and `stress_splicing_index < 0.4`
- `Unclassified` otherwise

Flags:
- `LOW_CONFIDENCE` if `confidence` is finite and `< 0.5`
- `LOW_SPLICE_SIGNAL` if `nnz < 50`
- `MISSING_METRICS` if any row metric is non-finite

Summary metrics in `summary.json`:
- Distribution stats for fidelity/stress: `median`, `p90`, `p99`
- Quantile rule: index `round((n-1)*q)` on sorted values
- `low_confidence_fraction = #cells(confidence < 0.5) / n_cells`
- `high_splice_noise_fraction = #cells(splice_junction_noise > 0.7) / n_cells`
- `splicing_instability` block:
  - panel metadata and thresholds
  - z-score reference medians/MAD
  - global `p50/p90` for `SOS`,`RLR`,`SII`
  - cluster stats (empty when cluster labels are unavailable)
  - missingness counters and panel coverage

## Caveats and Integration Note

- SOS/RLR/SII are transcriptional proxies, not direct molecular assays.
- RLR does not directly measure RNA:DNA hybrid occupancy.
- Interpretation is strongest when integrated with replication stress / DDR axes (e.g., `kira-nuclearqc`).
- Recommended downstream cross-axis predicate in `kira-organelle`:
  - `replication_stress_high && splicing_instability_high`

## Constants

- `EPS_STAGE3 = 1e-12` (entropy/dispersion stability in stage 3; the entropy z-score uses `EPS_ROBUST`)
- `EPS_ROBUST = 1e-6` (robust z-score denominator in stages 3, 4, 5, 8, 9, 10, 11)
- `EPS_STAGE11 = 1e-6` (splicing-noise denominator)
- `MAD_SCALE = 1.4826` (MAD to robust sigma scale)
- `SAT = 6.0` (stage 12 saturation bound)
- SIS weights: `0.35, 0.25, 0.25, 0.15`
- Missplicing burden weights: `0.35, 0.25, 0.25, 0.15`
- Confidence weights: `0.65` (fidelity), `0.35` (inverse penalty), penalty mix `0.4/0.3/0.3`
