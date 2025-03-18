process CORESNPFILTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/core-snp-filter:0.2.0--h4c94732_1' :
        'biocontainers/core-snp-filter:0.2.0--h3ab6199_2' }"

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("*.aln")      , emit: alignment
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    coresnpfilter \\
        --exclude_invariant \\
        ${alignment} \\
        ${args} > "${prefix}.aln"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coresnpfilter: \$( coresnpfilter --version 2>&1  | grep -o "v.*" )
    END_VERSIONS
    """


    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.aln"


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coresnpfilter: \$( coresnpfilter --version 2>&1  | grep -o "v.*" )
    END_VERSIONS

    """
    }
