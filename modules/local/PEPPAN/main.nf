process PEPPAN {
    tag "$meta.id"
    label 'process_medium'
    conda '$moduleDir/environment.yml'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cidgoh/peppan:latest' :
        'cidgoh/peppan:latest' }"

    input:
    tuple val(meta), path(gff)
    path reference_gff

    output:
    tuple val(meta), path("*.gff")   ,           emit: gff_peppan
    tuple val(meta), path("*.fna")   ,           emit: fna_peppan
    path "versions.yml"              ,           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ref_gff = reference_gff ? "-P $reference_gff" : ''
    def version = "1.0.5" // ------> update manually

    """
    PEPPAN -p $prefix \
        $ref_gff \
        $gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        peppan: $version
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.html
    touch ${prefix}.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        peppan: $version
    END_VERSIONS
    """
}
