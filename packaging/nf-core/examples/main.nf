#!/usr/bin/env nextflow
// Minimal workflow: one STARsolo Solo.out directory per sample from a
// samplesheet (`id,solo_out,metadata`), kira-spliceqc on each, MultiQC on
// the custom-content files.
//
//   nextflow run packaging/nf-core/examples/main.nf --samplesheet samples.csv --outdir results
nextflow.enable.dsl = 2

params.samplesheet = null
params.outdir = 'results'
params.reference = null
params.kira_spliceqc_subset = 'filtered'

include { KIRA_SPLICEQC } from '../modules/kira/spliceqc/main'

process MULTIQC {
    container 'multiqc/multiqc:v1.25'
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    input:
    path 'kira/*'
    output:
    path 'multiqc_report.html'
    script:
    """
    multiqc . --filename multiqc_report.html
    """
}

workflow {
    samples = Channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.id]
            def metadata = row.metadata ? file(row.metadata, checkIfExists: true) : []
            def reference = params.reference ? file(params.reference, checkIfExists: true) : []
            tuple(meta, file(row.solo_out, checkIfExists: true), metadata, reference)
        }
    KIRA_SPLICEQC(samples)
    KIRA_SPLICEQC.out.cells.map { meta, f -> f }.collectFile(storeDir: "${params.outdir}/cells")
    MULTIQC(KIRA_SPLICEQC.out.multiqc.map { meta, f -> f }.collect())
}
