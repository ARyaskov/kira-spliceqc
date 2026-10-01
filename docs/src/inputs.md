# Input levels

The tool uses whatever evidence the dataset carries and records the levels in
`summary.json.input.levels`.

| level | evidence | how it is found | enables |
| --- | --- | --- | --- |
| L0 | gene counts | `--input`: 10x MatrixMarket directory, `.h5ad`, `kira-organelle.bin` cache | expression signatures (`*_expr`), cell-cycle scores, cell QC flags |
| L1 | spliced / unspliced (/ ambiguous) counts | `spliced.mtx` + `unspliced.mtx` next to `matrix.mtx` (kb-python); STARsolo `Solo.out/Velocyto/<subset>`; `.h5ad` `layers/spliced`, `layers/unspliced`; `--layers PATH` | `unspliced_fraction`, `intron_retention_index`, `nuclear_fraction_flag`, `intron_retention_high` |
| L2 | junction counts | STARsolo `Solo.out/SJ/<subset>`; `--junctions DIR` | `cryptic_3ss_fraction`, `exon_skip_fraction`, `splice_site_shift`, `unannotated_junction_fraction` and their flags |

## Producing L1 and L2 with STARsolo

```bash
STAR --runThreadN 16 --genomeDir $INDEX \
  --readFilesIn R2.fastq.gz R1.fastq.gz --readFilesCommand zcat \
  --soloType CB_UMI_Simple --soloCBwhitelist $WHITELIST \
  --soloFeatures Gene Velocyto SJ --soloMultiMappers EM \
  --outSAMtype None --outFileNamePrefix sample/
```

`benchmarks/align_starsolo.sh` wraps this. kb-python (`kb count --workflow
lamanno` / `nac`) produces L1 as `spliced.mtx` / `unspliced.mtx`.

## Metadata

`--metadata PATH` (or `metadata.tsv[.gz]` next to the matrix, or `obs` of an
`.h5ad`) supplies cell annotations. Columns used:

- stratification: `--stratify-by COLUMN`, else the first `cell_type`-like,
  then `cluster`-like column;
- doublets: a boolean-like column (`predicted_doublet`, `doublet`,
  `is_doublet`, `scDblFinder.class`, `doublet_class`, `DF.classifications`).

## Cell QC

Cells below `--min-counts` (500) or `--min-genes` (200) are flagged
`LOW_DEPTH`; doublets are flagged `DOUBLET`. Both keep their raw metrics, get
no deviations or flags, and are excluded from every reference norm.
