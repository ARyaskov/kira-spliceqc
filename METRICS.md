# kira-spliceqc Metrics Specification

This document defines the metrics produced by `kira-spliceqc`, including formulas, constants, and classification rules.

Scope:
- standalone outputs: `spliceqc.tsv`, `spliceqc.json`
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

## Notation

Let:
- `count(g, c)` = raw UMI count for gene `g` in cell `c`.
- `libsize(c)` = sum of counts in cell `c` (clamped with `max(1, libsize)` where needed).
- `cp10k(g, c) = 1e4 * count(g, c) / max(1, libsize(c))`.
- `log_cp10k(g, c) = ln(1 + cp10k(g, c))`.
- `relu(x) = max(x, 0)`.
- `sigmoid(x) = 1 / (1 + exp(-x))`.

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
- `b_core = relu(-mean(z_u1, z_u2, z_sf3b))`
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
- if `MAD == 0`, z-score for finite values is `0`

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
- `cryptic_risk = sigmoid((x_sr_hnrnp + x_entropy + x_nmd) - 1.5)`

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

Row metrics in `kira-spliceqc/spliceqc.tsv`:
- `splice_fidelity_index = normalize01(sis)`
- `intron_retention_rate = normalize01(b_u12)`
- `exon_skipping_rate = normalize01(|axis_u2_u1|)`
- `alt_splice_burden = normalize01(missplicing_burden)`
- `splice_junction_noise = normalize01(burden_star)`
- `stress_splicing_index = normalize01(coupling_stress)` if stage 8 present, else `normalize01(imbalance)`

`normalize01(x)`:
- if finite and already in `[0,1]`, keep value
- otherwise apply logistic transform `1 / (1 + exp(-x))`, then clamp to `[0,1]`
- if non-finite, return `0`

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
- `LOW_CONFIDENCE` if `confidence < 0.5`
- `LOW_SPLICE_SIGNAL` if `nnz < 50`

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

- `EPS_STAGE3 = 1e-12` (entropy/dispersion stability in stage 3)
- `EPS_ROBUST = 1e-6` (robust z-score denominator in stages 4, 5, 8, 9, 10, 11)
- `EPS_STAGE11 = 1e-6` (splicing-noise denominator)
- `MAD_SCALE = 1.4826` (MAD to robust sigma scale)
- `SAT = 6.0` (stage 12 saturation bound)
- SIS weights: `0.35, 0.25, 0.25, 0.15`
- Missplicing burden weights: `0.35, 0.25, 0.25, 0.15`
- Confidence weights: `0.65` (fidelity), `0.35` (inverse penalty), penalty mix `0.4/0.3/0.3`
