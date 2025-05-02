process FAI2BED {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(fai)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path('versions.yml'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    """
    awk 'BEGIN {FS="\t"}; {print \$1 FS "0" FS \$2}' ${fai} > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$( awk --version | head -n1 | cut -f3 -d' ' | tr -d ',' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$( awk --version | head -n1 | cut -f3 -d' ' | tr -d ',' )
    END_VERSIONS
    """
}
