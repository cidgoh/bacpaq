process MEDAKA_VARIANT {
    tag "$meta.id"
    // label 'process_high'
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:2.0.1--py310he807b20_0' :
        'biocontainers/medaka:2.0.1--py310he807b20_0' }"

    input:
    tuple val(meta), path(input_reads)
    path(reference)

    output:
    tuple val(meta), path("*/${meta.id}*.vcf")    , emit: alignment_vcf
    tuple val(meta), path("*/*annotated.vcf")   , emit: annotated_vcf
    tuple val(meta), path("*/*sorted.vcf")      , emit: sorted_vcf
    tuple val(meta), path("*/*.bam.bai")        , emit: bam_bai
    tuple val(meta), path("*/*.bam")            , emit: bam
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    medaka_variant \\
        -i ${input_reads} \\
        -r ${reference} \\
        -o ${prefix} \\
        -t ${task.cpus} \\
        ${args} 

    mv ${prefix}/medaka.vcf ${prefix}/${prefix}.vcf
    mv ${prefix}/calls_to_ref.bam ${prefix}/${prefix}.bam
    mv ${prefix}/calls_to_ref.bam.bai ${prefix}/${prefix}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}/testing.vcf
    touch ${prefix}/testing_annotated.vcf
    touch ${prefix}/testing_sorted.vcf
    touch ${prefix}/testing.bam
    touch ${prefix}/testing.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS

    """
}
