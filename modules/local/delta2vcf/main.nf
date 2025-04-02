process DELTA2VCF {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mummer:3.23--pl5262h1b792b2_12' :
        'biocontainers/mummer:3.23--pl5262h1b792b2_12' }"

    input:
    tuple val(meta), path(delta)

    output:
    tuple val(meta), path("*.vcf")      , emit: vcf
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    delta2vcf.pl \\
        < ${delta} \\
        > "${prefix}.vcf"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$( perl --version 2>&1  | grep -o "v[0-9].*" | sed -E "s/).*//" )
    END_VERSIONS
    """


    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.vcf"


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$( perl --version 2>&1  | grep -o "v[0-9].*" | sed -E "s/\\).*//" )
    END_VERSIONS

    """
    }
