# kira-spliceqc (Python wrapper)

A thin wrapper around the `kira-spliceqc` binary for the scverse ecosystem:
it runs the tool and returns an `AnnData` whose `obs` carries every per-cell
column of `cells.tsv` and whose `uns["kira_spliceqc"]` carries `summary.json`.
The binary itself is not bundled: install it from bioconda
(`conda install -c bioconda kira-spliceqc`), from crates.io
(`cargo install kira-spliceqc`) or from a GitHub release, or point
`KIRA_SPLICEQC_BIN` at it.

```python
import kira_spliceqc as ksq

adata = ksq.run("data/pbmc3k/Solo.out/Gene/filtered", out="out/pbmc3k",
                metadata="data/pbmc3k/metadata.tsv")
adata.obs[["unspliced_fraction", "nuclear_fraction_flag", "cryptic_3ss_high"]].head()
adata.uns["kira_spliceqc"]["reference"]["mode"]
```

`run` accepts an `AnnData` too: it is written to a temporary `.h5ad`
(`layers["spliced"]` / `layers["unspliced"]` and `obs` are picked up as input
levels L1 and metadata) and the results are merged back into a copy.

`simulate`, `validate` and `reference_build` forward to the matching
subcommands and return the output paths; `kira-spliceqc-py ...` forwards any
command line to the binary.
