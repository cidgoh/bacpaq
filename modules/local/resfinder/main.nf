process RESFINDER {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::resfinder=4.1.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'staphb/resfinder:4.1.11' :
        'staphb/resfinder:4.1.11' }"
    
    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}/*"), optional: true, emit: output
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def set_acquired = acquired ? "--acquired" : ''
    def set_point = point ? "--point" : ''
    """
    run_resfinder.py \\
        -ifa ${fasta} \\
	-o ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        resfinder: \$( resfinder 2>&1 | grep '...:::' | sed 's/.*ResFinder v//;s/ .*//' )
    END_VERSIONS
    """
}
