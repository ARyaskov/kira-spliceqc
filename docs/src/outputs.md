# Outputs and pipeline integration

## Per-cell table (`cells.tsv`, `cells.json`)

One row per cell. Column families, in order:

- identity and QC: `cell_name`, `low_depth`, `doublet`;
- expression signatures (`*_expr`), present on any input;
- Tier A (L1): `spliced_umis`, `unspliced_umis`, `ambiguous_umis`,
  `unspliced_fraction` with `_ci_low` / `_ci_high` and `_dev`,
  `nuclear_fraction_flag`, `intron_retention_index`, `ir_gene_dispersion`,
  `ir_genes_used`, `intron_retention_index_dev`, `intron_retention_high`;
- Tier B (L2): `junction_umis`, `annotated_umis`,
  `unannotated_junction_fraction`, `cryptic_3ss_fraction` (+ `_dev`,
  `cryptic_3ss_high`), `exon_skip_fraction` (+ `_dev`, `exon_skip_high`),
  `splice_site_shift` (+ `_dev`, `splice_site_shift_high`), `site_groups_used`;
- cell cycle: `s_score_expr`, `g2m_score_expr`, `cell_cycle_phase`, `cycling`;
- experimental composites (`sis`, `class`, `SOS`, `RLR`, `SII`, flags), only
  with `--experimental-signatures` or in pipeline mode.

Undefined values are empty in TSV and `null` in JSON. Every `_dev` column is a
signed deviation in robust-z units; every `_high` / `_flag` column is boolean.

## Sample summary (`summary.json`)

`input` (levels, cells, species), `qc` (LOW_DEPTH / DOUBLET fractions),
`unspliced`, `intron_retention`, `junctions` (medians, flag fractions,
per-stratum references), `reference` (mode, strata, external file and the
metrics that used it), `cell_cycle`, the experimental `splicing_instability`
block, and `provenance` (command line, catalog and reference hashes, every
model constant, undefined-cell counts).

## Pipeline mode

`--run-mode pipeline` writes into `<OUT>/kira-spliceqc/` and adds:

- `spliceqc.tsv`: the pipeline-contract table (`barcode`, `regime`,
  `confidence`, `flags`, contract metrics); columns are versioned by
  `contract_version` in `pipeline_step.json`;
- `panels_report.tsv`: panel coverage per geneset;
- `pipeline_step.json`: manifest naming every artifact;
- `kira_spliceqc_mqc.json`: MultiQC custom content.

## MultiQC

```bash
multiqc ./out --filename splicing_qc.html
```

MultiQC finds `kira_spliceqc_mqc.json` under the search path and renders a
"kira-spliceqc" table with one row per sample: cells, LOW_DEPTH / DOUBLET /
cycling fractions, median unspliced fraction and damaged-cell fraction,
median intron retention index and IR-high fraction, Tier B medians and flag
fractions, and the reference mode. Tier columns appear only when the input
level was available.

## Nextflow / nf-core

`packaging/nf-core/modules/kira/spliceqc/` is an nf-core-style module (see
the [tutorial](tutorials/nf-core.md)).

## Python

`kira_spliceqc.run(...)` (package in `python/`) returns an `AnnData` with
`cells.tsv` merged into `obs` and `summary.json` in `uns["kira_spliceqc"]`.
