# Tutorial: finding SF3B1-like cells in CLL

Goal: in a chronic lymphocytic leukaemia sample, find cells whose 3' splice-site
usage matches the SF3B1 hotspot-mutant phenotype (cryptic acceptors 10-50 nt
upstream of canonical ones; Darman et al. 2015, Alsafadi et al. 2016), using
junction counts and an external reference from a healthy donor.

## Why an external reference

In a CLL sample most cells are tumour B cells. If the clone carries SF3B1
K700E, the B-cell stratum's own median cryptic usage *is* the mutant level and
a dataset-relative comparison flags nothing. The comparison has to be against
B cells of a control.

## 1. Align both samples with junction output

Run STARsolo with `--soloFeatures Gene Velocyto SJ` (see the PBMC tutorial)
for the CLL sample and for a healthy PBMC sample from the same chemistry.
Annotate both with the same labels (`B`, `T`, `NK`, `Mono`, ...).

## 2. Build the reference from the control

```bash
kira-spliceqc reference build \
  --input control/Solo.out/Gene/filtered --metadata control/metadata.tsv \
  --stratify-by cell_type --out ref_pbmc.json
```

## 3. Score the CLL sample against it

```bash
kira-spliceqc run \
  --input cll/Solo.out/Gene/filtered --metadata cll/metadata.tsv \
  --reference ref_pbmc.json --out out/cll --run-mode pipeline
```

`summary.json.reference.mode` is `external` and `external_metrics` lists
`unspliced_fraction`, `intron_retention_index` and `expression_signatures`.
Tier B deviations use the dataset's own strata (junction norms are not stored
in the reference yet), so for the cryptic fraction the key comparison is the
raw value and the between-stratum contrast, step 4.

## 4. Look at the cryptic fraction by cell type

```python
import pandas as pd
cells = pd.read_csv("out/cll/kira-spliceqc/cells.tsv", sep="\t", index_col="cell_name")
meta = pd.read_csv("cll/metadata.tsv", sep="\t", index_col="barcode")
cells = cells.join(meta)
cells.groupby("cell_type")["cryptic_3ss_fraction"].describe()
```

In an SF3B1-mutant CLL the B-cell median cryptic fraction is several-fold
above T cells and above control B cells (Nam et al. 2019 report the K700E
clone separating from wild-type cells on exactly these junctions). Compare
with the control's B cells:

```python
control = pd.read_csv("out/control/kira-spliceqc/cells.tsv", sep="\t", index_col="cell_name")
```

## 5. Per-cell calls within the tumour

For a subclonal mutation the tumour B cells themselves split: run the CLL
sample with `--stratify-by cell_type` (internal reference) and read
`cryptic_3ss_high` within B cells; the flagged fraction estimates the mutant
clone size. Confirm with genotyping (GoT amplicon, scDNA) before calling a
clone mutant: the tool reports a phenotype, not a genotype.

## 6. Corroborating signals

`intron_retention_index_dev` tends to rise in SF3B1 / SRSF2 / U2AF1 mutant
cells (Pellagatti et al. 2018); `exon_skip_fraction` changes in U2AF1 S34F.
Expression signatures (`sf3b`-related `*_expr`) are not evidence of the
mutation.
