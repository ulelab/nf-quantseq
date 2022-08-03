process RANDOM_PRIMING {
    tag "$meta.id"
    label 'process_medium'

    // conda (params.enable_conda ? 'bioconda::cutadapt=3.4' : null)
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/cutadapt:3.4--py39h38f01e4_1' :
    //     'quay.io/biocontainers/cutadapt:3.4--py37h73a75cf_1' }"
    container 'amchakrabarti/nf-quantseq'

    input:
    tuple val(meta), path(pos_bedgraph)
    tuple val(meta), path(neg_bedgraph)

    output:
    tuple val(meta), path("${prefix}.artefactsannotated.bed")       , emit: artefactsannotated_bed
    tuple val(meta), path("${prefix}.filteredunique.bed")           , emit: filteredunique_bed
    tuple val(meta), path('*.log')                                  , emit: log
    tuple val(meta), path('*.pdf')                                  , emit: plots
    path "versions.yml"                                             , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    remove_random_priming.R \
        --bg_pos $pos_bedgraph \
        --bg_neg $neg_bedgraph \
        --output ${prefix}.bed \
        --cores ${task.cpus} \
        $args \
        > remove_random_priming.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nf-quantseq: dev
    END_VERSIONS
    """
}
