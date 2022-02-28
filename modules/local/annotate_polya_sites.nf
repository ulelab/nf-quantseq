process ANNOTATE_POLYA_SITES {
    tag "$meta.id"
    label 'process_medium'

    // conda (params.enable_conda ? 'bioconda::cutadapt=3.4' : null)
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/cutadapt:3.4--py39h38f01e4_1' :
    //     'quay.io/biocontainers/cutadapt:3.4--py37h73a75cf_1' }"
    container 'amchakrabarti/nf-quantseq'

    input:
    tuple val(meta), path(bed)
    path gff3

    output:
    tuple val(meta), path('*.annotated.bed'), emit: bed
    tuple val(meta), path('*.log')          , emit: log
    path "versions.yml"                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    annotate_pas.R \
        --threads ${task.cpus} \
        --bed $bed \
        --gff3 $gff3 \
        > annotate_pas.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nf-quantseq: dev
    END_VERSIONS
    """
}