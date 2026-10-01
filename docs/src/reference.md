# Reference strata and external references

Every deviation (`*_dev`) and flag compares a cell with a reference: the
cells of its own stratum, within its library-size bin. Three modes exist and
`summary.json.reference.mode` records which ran:

| mode | when | norms come from |
| --- | --- | --- |
| `global` | no metadata column usable | all cells of the dataset |
| `stratified` | a cell-type / cluster column (or `--stratify-by`) | the cell's stratum; strata under 50 cells fold into `global` |
| `external` | `--reference ref.json` | a control dataset's strata, matched by the same metadata column |

Within a stratum, cells are split into up to 20 library-size bins of at least
50 cells each; the median and MAD of a bin define the robust z of its cells,
so the signatures are free of library-size correlation (enforced on the null
model). Proportions (unspliced fraction, cryptic usage, exon skipping) are
compared on the logit scale with a method-of-moments overdispersion so that
cells with few UMIs get wide uncertainty; the intron retention index carries a
per-cell standard error the same way.

A flag is set when the deviation is at least 3 in the expected direction and
its Benjamini-Hochberg adjusted p-value within the stratum is below 0.05.

## When the internal reference is wrong

A dataset-relative reference cannot see a shift that affects a whole stratum:
if every B cell of a sample carries an SF3B1 mutation, the B-cell stratum's
own median is the mutant level. Build a reference from a control dataset with
the same protocol and annotation and apply it:

```bash
kira-spliceqc reference build --input ./control/Solo.out/Gene/filtered \
  --metadata ./control/metadata.tsv --stratify-by cell_type --out ./ref.json
kira-spliceqc run --input ./sample/Solo.out/Gene/filtered \
  --metadata ./sample/metadata.tsv --reference ./ref.json --out ./out/sample
```

The file stores, per stratum, the unspliced-fraction and intron-retention
norms, pooled per-gene unspliced ratios, and depth-binned norms of every
geneset activity, the regulator entropy and the stage-15 cores.
`summary.json.reference.external_metrics` lists what used it. Cells whose
stratum is absent from the file fall back to `global`.

Protocol matters: nuclei have unspliced fractions around 0.5-0.7, whole cells
around 0.1-0.3. Never apply a whole-cell reference to nuclei or the reverse.
