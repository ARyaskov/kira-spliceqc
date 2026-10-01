# Interpretation rules

1. **A flag is a statistical call, not a diagnosis.** `*_high` /
   `nuclear_fraction_flag` mean: deviation of at least 3 robust-z from the
   cell's reference stratum and depth bin, with FDR 5 % within the stratum.
   On the null model at most 1 % of cells are flagged; in a real sample 1 %
   of flags are the expected false positives.
2. **Read the reference first.** `summary.json.reference.mode` tells whether
   a cell was compared with its own cell type (`stratified`), with the whole
   sample (`global`) or with a control dataset (`external`). A stratum-wide
   effect is invisible in the first two modes by construction.
3. **Measurements before proxies.** Tier A and Tier B columns measure
   splicing; `*_expr` columns measure the expression of splicing genes.
   A cell with high `cryptic_3ss_fraction_dev` and nothing else is a
   stronger finding than a cell with every `*_expr` column elevated.
4. **Damaged cells first.** `nuclear_fraction_flag` cells (and `LOW_DEPTH`,
   `DOUBLET`) should be removed or stratified before any biological reading;
   a damaged cell shifts every other metric.
5. **Check the support.** `ir_genes_used`, `site_groups_used`, `junction_umis`,
   `unspliced_fraction_ci_low/high`: a flag resting on few observations is
   wide by construction, but a cluster of such flags still deserves a look at
   the raw counts.
6. **Cell cycle is a covariate.** Expression signatures of the spliceosome
   rise in S / G2M; compare `cycling` fractions between groups before
   interpreting `spliceosome_core_expr`.
7. **Protocol decides the scale of the unspliced fraction.** Compare nuclei
   with nuclei and cells with cells; an external reference from another
   protocol is wrong, not conservative.
8. **Experimental composites are hypotheses.** `sis`, `SOS`, `RLR`, `SII` and
   the regimes exist for pipeline continuity; treat them as unvalidated until
   the tier 2-3 benchmarks are published.
9. **Pseudobulk before claiming a sample-level effect.** Per-cell flags
   aggregate to fractions (`summary.json`, the MultiQC table); a sample
   effect is a difference of fractions between samples, with replicates.
