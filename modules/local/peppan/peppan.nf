process PEPPAN_PEPPAN {
    tag "$meta.id"
    label 'process_medium'
    conda 'assets/peppan/environment.yml'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cidgoh/peppan:latest' :
        'cidgoh/peppan:latest' }"

    input:
    tuple val(meta), path(gff)
    path (reference_gff)

    output:
    tuple val(meta), path("*.gff"), path("*.fna")   ,             emit: gff_peppan
    path "versions.yml"                             ,             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def version = "1.0.5" // ------> update manually
    """
    PEPPAN -p $prefix \
        -P $reference_gff \
        *gff.gz
    
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

params.gff_reference = null

workflow {
    
    ch_gff              =   Channel.fromPath(params.gff_peppan)
    ch_reference_gff    =   Channel.fromPath(params.gff_reference)
    ch_gff
        .collect()
        .map {tuple([ id:'test' ], it) }
        .set{ch_all}

     PEPPAN_PEPPAN(ch_all,
            ch_reference_gff)
}