# Tutorial: running inside nf-core/scrnaseq

nf-core/scrnaseq with `--aligner star` produces `Solo.out` directories per
sample; `kira-spliceqc` adds a splicing QC step and a MultiQC section.

## Option A: the standalone example workflow

`packaging/nf-core/examples/main.nf` runs the module on a samplesheet and
collects the MultiQC report:

```csv
id,solo_out,metadata
patient1,results/star/patient1/patient1.Solo.out,meta/patient1.tsv
patient2,results/star/patient2/patient2.Solo.out,
```

```bash
nextflow run packaging/nf-core/examples/main.nf \
  --samplesheet samples.csv --outdir results/spliceqc \
  --kira_spliceqc_subset filtered -profile docker
```

Set `--reference ref.json` to apply an external reference to every sample.

## Option B: the module in your own pipeline

Copy `packaging/nf-core/modules/kira/spliceqc/` into `modules/local/` (or,
once it is merged into nf-core/modules, `nf-core modules install kira/spliceqc`)
and include it:

```groovy
include { KIRA_SPLICEQC } from '../modules/local/kira/spliceqc/main'

workflow {
    ch_solo = STARSOLO.out.solo_out.map { meta, dir -> tuple(meta, dir, [], []) }
    KIRA_SPLICEQC(ch_solo)
    ch_multiqc_files = ch_multiqc_files.mix(KIRA_SPLICEQC.out.multiqc.map { it[1] })
    ch_versions = ch_versions.mix(KIRA_SPLICEQC.out.versions)
}
```

The process takes `tuple(meta, input, metadata, reference)`; pass `[]` for an
absent metadata table or reference. `task.ext.args` forwards extra flags
(`--stratify-by cluster`, `--min-counts 1000`, `--extended`), `task.ext.prefix`
sets the output directory name.

## Resources

One sample of 10 000 cells and 30 000 genes runs in well under a minute on
4 threads with a few hundred MB of memory; `label 'process_medium'` is
generous. The container is `biocontainers/kira-spliceqc`; the conda
environment is `bioconda::kira-spliceqc`.
