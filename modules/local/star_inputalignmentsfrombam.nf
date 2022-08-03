process STAR_INPUTALIGNMENTSFROMBAM {
    tag "$meta.id"
    label 'process_high'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda (params.enable_conda ? 'bioconda::star=2.7.9a' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.9a--h9ee0642_0' :
        'quay.io/biocontainers/star:2.7.9a--h9ee0642_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("$prefix*.Unique.str1*")        , emit: unique_coverage_pos
    tuple val(meta), path("$prefix*.Unique.str2*")        , emit: unique_coverage_neg
    tuple val(meta), path("$prefix*.UniqueMultiple.str1*"), emit: unique_multiple_coverage_pos
    tuple val(meta), path("$prefix*.UniqueMultiple.str2*"), emit: unique_multiple_coverage_neg
    path "versions.yml"                                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    STAR \\
        --runMode inputAlignmentsFromBAM \\
        --inputBAMfile $bam \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}
