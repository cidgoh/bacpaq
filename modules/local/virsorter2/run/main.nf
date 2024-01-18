process VIRSORTER2_RUN {
    tag "$meta.id"
    label 'process_medium'
    conda "bioconda::virsorter=2.2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://jiarong/virsorter:2.2.3' :
        'jiarong/virsorter:2.2.3' }"

    input:
    tuple val(meta), path(fasta)
    path(db)

    output:
    path "final-viral-boundary.tsv" , emit: boundary                                     
    path "final-viral-combined.fa"  , emit: fasta
    path "final-viral-score.tsv"    , emit: score
    path "log"                      , emit: log
    path "config.yaml"              , emit: config
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def db_path  = db ?: '/db'
    def prefix   = task.ext.prefix ?: meta.id
    """
    virsorter run \
        -d ${db_path} \
        -j ${task.cpus} \
        -w . \
        -i ${fasta} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            virsorter --version | tail -n1 | sed 's/.*version //g')
    END_VERSIONS
    """
}