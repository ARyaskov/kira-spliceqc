// nf-core-style module for kira-spliceqc. Input: a STARsolo `Solo.out`
// directory (Gene + optional Velocyto / SJ subdirectories) or any 10x
// MatrixMarket directory / .h5ad file; optional metadata table and external
// reference. Output: the pipeline-contract directory with cells.tsv,
// summary.json, spliceqc.tsv and the MultiQC custom-content file.
process KIRA_SPLICEQC {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kira-spliceqc:0.5.0--h9ee0642_0' :
        'biocontainers/kira-spliceqc:0.5.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(input), path(metadata), path(reference)

    output:
    tuple val(meta), path("${prefix}/kira-spliceqc/cells.tsv")      , emit: cells
    tuple val(meta), path("${prefix}/kira-spliceqc/cells.json")     , emit: cells_json, optional: true
    tuple val(meta), path("${prefix}/kira-spliceqc/summary.json")   , emit: summary
    tuple val(meta), path("${prefix}/kira-spliceqc/spliceqc.tsv")   , emit: contract
    tuple val(meta), path("${prefix}/kira-spliceqc/*_mqc.json")     , emit: multiqc
    tuple val(meta), path("${prefix}/kira-spliceqc/pipeline_step.json"), emit: step
    path "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def meta_arg = metadata ? "--metadata ${metadata}" : ''
    def ref_arg = reference ? "--reference ${reference}" : ''
    // A STARsolo Solo.out directory: point at Gene/<subset>; Velocyto and SJ
    // siblings are detected automatically.
    def input_arg = input.isDirectory() && input.resolve('Gene').exists() ?
        "${input}/Gene/${params.kira_spliceqc_subset ?: 'filtered'}" : "${input}"
    """
    kira-spliceqc run \\
        --input ${input_arg} \\
        --out ${prefix} \\
        --run-mode pipeline \\
        --threads ${task.cpus} \\
        ${meta_arg} \\
        ${ref_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kira-spliceqc: \$(kira-spliceqc --version | sed 's/kira-spliceqc //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/kira-spliceqc
    touch ${prefix}/kira-spliceqc/cells.tsv
    echo '{}' > ${prefix}/kira-spliceqc/summary.json
    touch ${prefix}/kira-spliceqc/spliceqc.tsv
    echo '{}' > ${prefix}/kira-spliceqc/kira_spliceqc_mqc.json
    echo '{}' > ${prefix}/kira-spliceqc/pipeline_step.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kira-spliceqc: \$(kira-spliceqc --version | sed 's/kira-spliceqc //')
    END_VERSIONS
    """
}
