process AWK {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::gawk=5.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'quay.io/biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(input)
    val command
    val post_command

    output:
    tuple val(meta), path("${prefix}.${extension}"), emit: file
    path  "versions.yml"                           , emit: versions

    script:
    prefix    = task.ext.prefix ?: "${meta.id}"
    extension = task.ext.extension ?: 'txt'
    def args  = task.ext.args ?: ''
    """
    awk $args '$command' $input $post_command > ${prefix}.${extension}
    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        awk: \$(awk -W version | head -n 1 | egrep -o "([0-9]{1,}\\.)+[0-9]{1,}")
    END_VERSIONS
    """
}
