process VIRSORTER2_RUN {
    tag "$meta.id"
    label 'process_medium'
    conda "bioconda::virsorter=2.2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://jiarong/virsorter:2.2.3' :
        'jiarong/virsorter:2.2.3' }"

    input:
    tuple val(meta), path(fasta)

    output:
    path "${prefix}/final-viral-boundary.tsv" , emit: boundary                                     
    path "${prefix}/final-viral-combined.fa"  , emit: fasta
    path "${prefix}/final-viral-score.tsv"    , emit: score
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def db       = params.virsorter_db ?: '/db'
    def prefix       = task.ext.prefix ?: '.'
    """
    virsorter run \
        -d ${db} \
        -j ${task.cpus} \
        -w ${prefix} \
        -i ${fasta} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            virsorter --version | tail -n1 | sed 's/.*version //g')
    END_VERSIONS
    """
}