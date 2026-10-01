# Tutorial: QC of a PBMC run from STARsolo

Goal: add splicing QC to a 10x PBMC sample, remove damaged cells, and look at
the result in MultiQC. Any 10x 3' or 5' library works; the example uses the
10x "pbmc_1k_v3" FASTQs.

## 1. Align with STARsolo, all three feature types

```bash
STAR --runThreadN 16 --genomeDir ./index/GRCh38 \
  --readFilesIn pbmc_1k_v3_S1_L001_R2_001.fastq.gz pbmc_1k_v3_S1_L001_R1_001.fastq.gz \
  --readFilesCommand zcat \
  --soloType CB_UMI_Simple --soloCBwhitelist 3M-february-2018.txt \
  --soloUMIlen 12 --soloFeatures Gene Velocyto SJ --soloMultiMappers EM \
  --soloCellFilter EmptyDrops_CR --outSAMtype None --outFileNamePrefix pbmc/
```

`pbmc/Solo.out/` now holds `Gene/filtered`, `Velocyto/filtered` (spliced /
unspliced / ambiguous) and `SJ/raw` (junction counts per barcode).

## 2. Annotate cell types (any method)

A metadata table with `barcode` and `cell_type` columns lets the tool compare
every cell with its own type. From Scanpy, after clustering and marker-based
labelling:

```python
adata.obs[["cell_type", "predicted_doublet"]].rename_axis("barcode").to_csv("pbmc/metadata.tsv", sep="\t")
```

## 3. Run

```bash
kira-spliceqc run \
  --input pbmc/Solo.out/Gene/filtered \
  --metadata pbmc/metadata.tsv \
  --out out/pbmc --run-mode pipeline
```

The log reports the input levels (`L0, L1, L2`), the strata (cell types with
\\(\ge 50\\) cells; smaller ones fold into `global`) and the stage timings.

## 4. Read the summary

```bash
jq '.input.levels, .qc, .unspliced | {median, nuclear_fraction_flag_fraction}' out/pbmc/kira-spliceqc/summary.json
```

Expect a median unspliced fraction around 0.15-0.25 for whole PBMCs and a
damaged-cell fraction of a few per cent at most. A damaged fraction above
10 % points at a stressed sample (long processing time, freeze-thaw).

## 5. Filter and stratify in Python

```python
import pandas as pd
cells = pd.read_csv("out/pbmc/kira-spliceqc/cells.tsv", sep="\t", index_col="cell_name")
keep = ~(cells.low_depth | cells.doublet | cells.nuclear_fraction_flag)
adata = adata[keep.reindex(adata.obs_names).fillna(False).values].copy()
adata.obs = adata.obs.join(cells[["unspliced_fraction", "intron_retention_index", "cryptic_3ss_fraction", "cycling"]])
```

(`kira_spliceqc.run(...)` from the Python package does the join for you.)

## 6. MultiQC

```bash
multiqc out/ --filename out/splicing_qc.html
```

The "kira-spliceqc" table shows one row per sample; with several samples the
damaged-cell and IR-high fractions become the first thing to compare.
