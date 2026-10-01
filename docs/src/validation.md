# Validation and benchmarks

Three tiers; the first runs in CI on every commit.

## Tier 1: simulation with known truth

`kira-spliceqc simulate` writes a Poisson dataset (two cell types, filler
genes, panel genes, spliced / unspliced layers, a STARsolo-style junction
matrix) and spikes four effects into disjoint cell subsets:

| effect | truth column | metric | flag |
| --- | --- | --- | --- |
| cryptic 3' splice-site usage | `truth_cryptic` | `cryptic_3ss_fraction_dev` | `cryptic_3ss_high` |
| global intron retention | `truth_ir` | `intron_retention_index_dev` | `intron_retention_high` |
| damaged cells (low nuclear fraction) | `truth_damaged` | `unspliced_fraction_dev` (negative) | `nuclear_fraction_flag` |
| exon skipping | `truth_skip` | `exon_skip_fraction_dev` | `exon_skip_high` |

`tests/validation_tier1.rs` requires AUROC \\(\ge 0.95\\), flag FPR
\\(\le 1\\,\%\\) and recall \\(\ge 0.9\\) per effect; the current result is
AUROC 1.0 and FPR below 0.5 %. The zero-effect case is the null model
(`tests/null_model.rs`): production flags \\(\le 1\\,\%\\) of cells and
\\(|\rho_{Spearman}(metric, library\ size)| \le 0.1\\) for every metric.

```bash
./benchmarks/run_tier1.sh ./bench-out 2000 24301
```

## Tier 2: positive controls on public data

Specified in `benchmarks/datasets.tsv` with acceptance criteria: GoT CLL and
MDS cohorts with per-cell SF3B1 / SRSF2 / U2AF1 genotypes (GSE117063),
genome-wide Perturb-seq with splicing-factor knockdowns (GSE202088), SF3b
inhibitor dose series, matched nuclei and whole cells. **Not executed yet**:
they need STAR alignments of the public FASTQs.

## Tier 3: agreement with published tools

velocyto / scVelo (unspliced fraction, Spearman \\(\ge 0.99\\)), DropletQC
(damaged-cell calls, kappa \\(\ge 0.7\\)), SpliZ (splice-site shift, per-gene
Spearman \\(\ge 0.8\\)), IRFinder (pseudobulk IR, Spearman \\(\ge 0.7\\)).
**Not executed yet.**

Processed matrices, results and reproduction scripts will be published on
Zenodo with a DOI once tiers 2-3 run; other tools are invited to add
themselves to the comparison (open a "New validation dataset" issue).
