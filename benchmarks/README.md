# kira-spliceqc benchmark package

Three validation tiers, each reproducible from this directory. Tier 1 runs
anywhere in minutes; tiers 2 and 3 need STAR, the reference genome and the
public FASTQs listed in `datasets.tsv`.

## Tier 1 — simulations with known truth (runs in CI)

```bash
./benchmarks/run_tier1.sh ./bench-out
```

`kira-spliceqc simulate` writes a Poisson dataset with two cell types, spliced /
unspliced layers and a STARsolo-style junction matrix, then spikes four effects
into disjoint cell subsets (cryptic 3' splice-site usage, global intron
retention, damaged cells, exon skipping). `kira-spliceqc run` scores it and
`kira-spliceqc validate` reports AUROC / AUPRC per metric and precision /
recall / FPR per flag, overall and per cell type. The null model
(`tests/null_model.rs`) is the zero-effect special case.

## Tier 2 — positive controls on public data

| dataset | accession | truth | metric under test | expected |
| --- | --- | --- | --- | --- |
| CLL, Genotyping of Transcriptomes (Nam et al. 2019 Nature) | GSE117063 | SF3B1 K700E per cell (GoT amplicon) | `cryptic_3ss_fraction_dev`, `cryptic_3ss_high` | AUROC >= 0.75 mutant vs wild-type within patient |
| MDS CD34+ with SF3B1 / SRSF2 / U2AF1 (GoT, Nam et al. 2019) | GSE117063 | genotype per cell | `cryptic_3ss_fraction_dev`, `intron_retention_index_dev` | effect sign and order as in bulk (Pellagatti et al. 2018 Blood) |
| K562 Perturb-seq genome-wide (Replogle et al. 2022 Cell) | GSE202088 | CRISPRi target per cell | `intron_retention_index_dev`, `_expr` signatures | snRNP / SF3B knockdowns above controls, AUROC >= 0.7 |
| SF3b inhibitor dose series (pladienolide B / E7107 / H3B-8800, cell lines) | see `datasets.tsv` | dose and time | `intron_retention_index`, `exon_skip_fraction` | monotone dose response, Spearman >= 0.8 across samples |
| Mixed nuclei / whole cells of one tissue | see `datasets.tsv` | protocol | `unspliced_fraction` | AUROC >= 0.95 protocol separation; none after `--protocol` stratification |

Workflow per dataset:

```bash
./benchmarks/align_starsolo.sh <fastq_dir> <star_index> <whitelist> <out>   # Gene + Velocyto + SJ
./benchmarks/run_dataset.sh <out>/Solo.out/Gene/raw <truth.tsv> <bench-out>  # run + validate
```

`align_starsolo.sh` requests `--soloFeatures Gene Velocyto SJ`, which gives all
three input levels (L0, L1, L2) in one pass. Truth tables are one row per
barcode with boolean columns (e.g. `truth_sf3b1`) and an optional `cell_type`;
pass the scoring pairs to `validate` with `--pair truth_sf3b1:cryptic_3ss_fraction_dev:cryptic_3ss_high`.

## Tier 3 — agreement with published tools

| ours | reference tool | criterion |
| --- | --- | --- |
| `unspliced_fraction` | velocyto / scVelo (same layers) | Spearman >= 0.99 per cell |
| `nuclear_fraction_flag` | DropletQC | Cohen's kappa >= 0.7 |
| `splice_site_shift` | SpliZ on the same BAM | Spearman >= 0.8 per gene |
| pseudobulk `intron_retention_index` | IRFinder on the same reads | Spearman >= 0.7 per gene |

## Status

Tier 1 passes (see `tests/validation_tier1.rs` and the CI log). Tiers 2 and 3
are specified here and in `datasets.tsv` but have not been executed yet: they
need STAR alignments of the public FASTQs, which this repository does not ship.
Processed matrices and results will be published on Zenodo with a DOI once run.
