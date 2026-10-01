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
| `external` | `--reference ref.json` (built by `kira-spliceqc reference build` on a control dataset) | cells are assigned to the reference strata by the reference's metadata column (unmatched -> `global`); Tier A deviations (`unspliced_fraction_dev`, `intron_retention_index_dev`) and the per-gene intron-retention ratios use the file's norms; expression signatures (geneset activity, regulator entropy, stage-15 cores) are standardized against the file's depth-binned norms when it has them (`expression_signatures` in `summary.json.reference.external_metrics`), otherwise they stay dataset-relative |

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

Reference file (`ref.json`, `format = kira-spliceqc-reference`, `version = 1`):
per stratum `unspliced_fraction {median_logit, tau2, median, n_defined}`,
`intron_retention_index {median, tau2, n_defined}` and
`gene_unspliced_ratio {symbol: p_gs}` and `expression {name: {edges, norms}}`,
the per-depth-bin median / MAD (`norms[i] = {median, mad, n}`; a cell takes
the first bin whose `edges[i]` is at or above its library size) of every
catalog geneset's raw activity, of the regulator entropy (`regulator_entropy`)
and of the stage-15 cores (`splice_core`, `rbp_core`, `rloop_resolve_core`,
`conflict_risk_core`, `nmd_core`), keyed by the same names the internal
standardization uses; the first stratum is always `global`. A stratum with no
defined cells, or a depth bin whose value is undefined or constant, is left
out rather than stored as NaN; cells landing there get NaN.
With an external reference the target run's own norms are still computed and
reported (`summary.json.unspliced.strata`), while `_dev` values and flags use the
file's norms, so a whole stratum that shifted relative to the control shows up
(it is invisible to a dataset-relative reference by construction).

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

## Tier B: Junction Metrics (Stage 19, requires input level L2)

Annotation is derived from the junction matrix itself: every junction the
aligner marked `annotated` defines an annotated donor and acceptor (strand-aware,
STAR coordinates = 1-based intron bounds). Junctions with unknown strand are
counted but never classified.

Classification of a junction `D -> A`:
- **cryptic 3' splice site**: unannotated, `D` is an annotated donor, and an annotated
  acceptor `A*` of the same donor lies `CRYPTIC_MIN..=CRYPTIC_MAX = 10..=50` nt
  downstream of `A` in transcript direction (i.e. `A` is 10-50 nt upstream of `A*`).
  The canonical partner is `D -> A*`. SF3B1 hotspot mutants (K700E, K666N, ...) shift
  branch-point recognition and select exactly such acceptors
  (Darman et al. 2015 Cell Reports; Alsafadi et al. 2016 Nature Communications).
- **exon skipping**: annotated `D -> A1` and `D2 -> A` exist with `A1` before `D2`
  (the junction skips at least one annotated exon; annotated skip junctions count).
  Inclusion partners are all such `D -> A1` and `D2 -> A` junctions
  (rMATS skipped-exon event definition, Shen et al. 2014 PNAS).
- **novel** otherwise.

Per cell `c` (undefined when `junction_umis < MIN_JUNCTION_UMIS = 200`):
- `junction_umis`, `annotated_umis`
- `unannotated_junction_fraction = (junction_umis - annotated_umis) / junction_umis`
- `cryptic_3ss_fraction = cryptic / (cryptic + canonical_partner)`; undefined when the
  denominator is below `MIN_RATIO_UMIS = 20`
- `exon_skip_fraction = skip / (skip + inclusion_partner / 2)`; same floor. An
  included exon is supported by two junctions and a skipped one by a single
  junction, so inclusion UMIs are halved (the complement of the rMATS
  junction-count PSI)
- `splice_site_shift` (after SpliZ, Olivieri et al. 2022 Nature Methods): for every
  donor with >= 2 acceptors and every acceptor with >= 2 donors, each UMI carries the
  rank of its partner site in transcript direction; with `r_c` the cell's mean rank at
  the site (>= `MIN_SITE_UMIS = 3` UMIs) and `r_s`, `v_s` the per-UMI mean and variance
  of the reference stratum, `z = (r_c - r_s) / sqrt(v_s / n_c)`; the score is
  `median_sites |z| / 0.6745` over >= `MIN_SITE_GROUPS = 5` sites (~1 under the null)
- `*_dev` and `*_high`: logit deviations with overdispersion (fractions) or robust z
  within stratum and junction-depth bin (shift, whose raw score rises with junction
  depth), flags at `dev >= 3` with BH-adjusted p < 0.05

Caveats: 3' 10x libraries cover few junctions per cell, so most cells may fall
below 200 junction UMIs; aggregate by cluster for such data (planned). Cryptic
detection needs the canonical acceptor to be used in the same dataset.

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
- `r_gc = log2(IR_gc / p_gs)` (log2 units; 0 = at the stratum reference, +1 = twice the
  reference unspliced ratio)
- `Var(r_gc) = n (1 - p_gs) / ((n + K)^2 p_gs ln^2 2)` (delta method), weights
  `w_gc = min(1 / Var(r_gc), 1 / Var at n = WEIGHT_CAP_UMIS = 50)` so no single gene dominates
- `intron_retention_index = sum_g w_gc r_gc / sum_g w_gc`; `se_c = 1 / sqrt(sum_g w_gc)`
- `ir_gene_dispersion = MAD_{g in G_c} r_gc`: small when every gene shifts together
  (global retention), large for gene-specific changes
- `ir_genes_used = |G_c|`
- `intron_retention_index_dev = (index_c - median_bin) / sqrt(se_c^2 + tau2_bin)` within the
  cell's stratum and layer-depth bin (`tau2` = method-of-moments overdispersion of the bin);
  `intron_retention_high` = `dev >= 3` and BH-adjusted two-sided normal p < 0.05

The raw index keeps a small depth bias (the log of a shrunk ratio is biased low
at few UMIs per gene); compare cells through `_dev`, use the raw index only for
effect sizes.

Adapted from the IRFinder intron-retention ratio (Middleton et al. 2017 Genome
Biology) to sparse per-cell counts. 3' libraries confound unspliced signal with
internal priming; the per-gene reference absorbs gene-specific background, so
the index measures departure from the stratum, not absolute retention.

## Geneset Activity (Stage 2)

Raw panel mean for geneset `S` and cell `c`:
- `A_raw(S, c) = mean_{g in S}( ln(1 + 1e4 * count(g, c) / max(1, libsize(c))) )`

Control-gene background (Tirosh et al. 2016 Science; Seurat `AddModuleScore`):
- gene mean `m_g = mean_c ln(1 + cp10k(g, c))` over all cells
- control pool = genes not in any splicing panel (catalog genesets and the stage-15
  panels) with `m_g > 0`, ranked by `m_g`
- for every panel gene the `CONTROLS_PER_GENE = 50` pool genes nearest in `m_g`;
  the panel's control set `C(S)` is their union (a panel gene is never its own control)
- `A(S, c) = A_raw(S, c) - mean_{g in C(S)} ln(1 + cp10k(g, c))`

A panel with an empty control set (every gene excluded) falls back to the raw
mean with a warning.

Standardization (`standardize_activity`): every panel score is turned into a
robust z-score against the cells of the same reference stratum **and**
library-size bin (`robust_z_by_stratum_and_depth`, up to 20 bins of at least
50 cells per stratum). The mean *and* the variance of a sparse panel score
depend on depth, so a stratum-wide reference removes only the mean. A bin with
a zero MAD leaves its cells undefined (no silent zeros). Stages 4, 5, 8, 9, 10
and 13 consume the standardized matrix; stage 11 (noise) uses the raw one;
stage 3 (`z_entropy`) and stage 15 (panel cores) apply the same
control-gene subtraction and depth-binned standardization.

On a Poisson null model this brings |Spearman(metric, libsize)| from 0.4-0.7
down to below 0.1 for `spliceosome_core_expr`, `SOS`,
`spliceosome_imbalance_expr` and `missplicing_burden_expr`
(`tests/null_model.rs`).

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
- `SPLICEQC_INSTABILITY_PANEL_V2`

Panels (human symbols, stable order):
- Core spliceosome / snRNP load proxy:
  - `SNRPB,SNRPD1,SNRPD2,SNRPD3,SNRPE,SNRPF,SNRPG,SF3A1,SF3A2,SF3A3,SF3B1,SF3B2,SF3B3,SF3B4,SF3B5,PRPF3,PRPF4,PRPF6,PRPF8,PRPF19,U2AF1,U2AF2`
- Splicing regulation / stress-sensitive RBPs:
  - `HNRNPA1,HNRNPA2B1,HNRNPC,HNRNPK,SRSF1,SRSF2,SRSF3,SRSF6,SRSF7,RBM39,RBM10,RBM17`
- R-loop resolution (protective axis):
  - `SETX,DDX5,DDX21,DHX9,RNASEH1,RNASEH2A,RNASEH2B,RNASEH2C,BRCA1,BRCA2`
- Transcription-replication conflict risk (optional):
  - `TOP1,TOP2B,POLR2A,SUPT5H,SUPT6H` (`TOP2A` removed in v0.5, panel V2: G2/M marker)
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

Flags (recalibrated in v0.5; the fixed cut-offs SOS >= 2.0 / RLR >= 1.5 /
SII >= 2.0 fired on 4 % of null cells):
- signed composites `SOS`, `RLR_signed = 0.7 * (-Z(rloop)) + 0.3 * Z(conflict)` (or
  `-Z(rloop)`), `SII_signed = 0.6 * SOS + 0.4 * (-Z(nmd))` (or `SOS`) are standardized
  within the reference stratum (`sos_dev`, `rlr_dev`, `sii_dev`)
- `splice_overload_high`: `sos_dev >= 3` and BH-adjusted p < 0.05 within the stratum
- `rloop_risk_high`: same rule on `rlr_dev`
- `splicing_instability_high`: same rule on `sii_dev`
- `genome_instability_splicing_flag`: `splice_overload_high && rloop_risk_high`
- the reported `SOS`/`RLR`/`SII` values keep the relu-based definitions above

NaN behavior:
- flags evaluate to `false` for non-finite scores
- missingness counters are emitted in JSON outputs

Optional junction-aware mode (future placeholder):
- `JE` (junction entropy)
- `IRB` (intron retention burden)
- `CSP` (cryptic splicing proxy)

## Cell-Cycle Scores (Stage 18, confounder annotation)

Gene lists: Tirosh et al. 2016 Science S-phase (43 genes) and G2/M (54 genes)
as in Seurat `cc.genes.updated.2019`; legacy symbols `MLF1IP`, `RPA2`,
`FAM64A`, `HN1` are recognised.

- `s_score_expr = mean_{g in S}(log1p cp10k) - mean_{g in ctrl(S)}(log1p cp10k)` with the
  control pool of "Geneset Activity" (cell-cycle genes are excluded from every control set);
  `g2m_score_expr` likewise
- undefined when fewer than 10 genes of either list map to the matrix
- `cell_cycle_phase` (Seurat `CellCycleScoring` rule): `S` if `s > g2m` and `s > 0`, `G2M`
  if `g2m > s` and `g2m > 0`, else `G1`; empty when undefined
- `cycling = phase in {S, G2M}`; the pipeline contract adds the `CYCLING` flag

Interpretation: spliceosome, TOP2A/BRCA1/BRCA2 and other panel genes rise in
S/G2M, so a cycling cell's expression signatures may reflect proliferation
rather than splicing regulation. The Seurat rule over-calls S/G2M on
non-cycling tissue (any positive noise counts), so treat the phase as an
annotation, not a QC verdict; `summary.json.cell_cycle` reports the fractions.

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
- `CYCLING` if `cell_cycle_phase` is `S` or `G2M`
- `LOW_DEPTH` if `libsize < min_counts` (default 500) or `nnz < min_genes` (default 200)
- `DOUBLET` if a doublet metadata column (`predicted_doublet`, `doublet`, `is_doublet`,
  `scDblFinder.class`, `doublet_class`, `DF.classifications`) is truthy (`true`, `1`, `yes`, `doublet`)

`LOW_DEPTH` and `DOUBLET` cells are excluded from every reference norm and get
undefined `_dev` values and flags (`cells.tsv` columns `low_depth`, `doublet`;
`summary.json.qc.low_depth_fraction` / `doublet_fraction`).

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
