# What the tool does not do

- **No alignment, no annotation.** It reads count matrices. Junction
  annotation is inferred from the annotated junctions STARsolo reports; a
  GTF is never parsed, so "annotated" means "annotated in the STAR index".
- **No isoform quantification or differential splicing.** No PSI per event,
  no gene-level tests, no splice-graph inference; SpliZ, scQuint, rMATS,
  MAJIQ or LeafCutter remain the tools for those questions.
- **No per-gene intron retention calls per cell.** The intron retention index
  is a per-cell summary over genes; per-gene IR needs pseudobulks and a tool
  such as IRFinder.
- **No mutation calling.** A cell with high cryptic 3' splice-site usage is
  SF3B1-*like*; genotype comes from GoT, scDNA-seq or amplicon data.
- **No RNA velocity.** Spliced / unspliced layers are used for QC only.
- **No doublet or ambient detection.** Doublet calls are read from metadata;
  ambient RNA is not modelled.
- **No batch correction.** Strata are compared with themselves; batch is a
  stratification choice (`--stratify-by`) or an external reference, not a
  correction.
- **No machine-learned scores.** Every number is a closed-form statistic with
  constants listed in `summary.json.provenance`; this is a feature.
- **No sample mode yet.** `--mode sample` returns an error; aggregate numbers
  come from `summary.json` and the MultiQC table.
- **Not validated on real positive controls yet.** Tier 1 (simulation) is
  enforced in CI; tiers 2-3 are specified in `benchmarks/` and pending.
