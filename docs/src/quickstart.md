# Quick start

A STARsolo run with `--soloFeatures Gene Velocyto SJ` gives all three input
levels; point the tool at the gene matrix and the siblings are detected:

```bash
kira-spliceqc run \
  --input ./sample/Solo.out/Gene/filtered \
  --metadata ./sample/metadata.tsv \
  --out ./out/sample \
  --run-mode pipeline
```

Outputs land in `./out/sample/kira-spliceqc/`:

| file | content |
| --- | --- |
| `cells.tsv` | one row per cell: raw metrics, `_dev` deviations, `_high` / `_flag` calls |
| `summary.json` | sample-level numbers, reference strata, provenance |
| `spliceqc.tsv` | pipeline-contract table (regime, confidence, flags) |
| `kira_spliceqc_mqc.json` | MultiQC custom content (`multiqc ./out/sample`) |

Without layers or junctions the tool still runs on gene counts alone and
reports the expression signatures (`*_expr`), which are proxies, not
measurements (see [Interpretation rules](interpretation.md)).

Check the installation on a simulated dataset with known truth:

```bash
kira-spliceqc simulate --out ./sim
kira-spliceqc run --input ./sim --junctions ./sim/sj --out ./sim-out --run-mode pipeline
kira-spliceqc validate --run ./sim-out --truth ./sim/truth.tsv --out ./validation.json
```

`validation.md` next to `validation.json` lists AUROC per metric and
precision / recall / FPR per flag; every spiked effect should be at AUROC 1.0.
