process VIRSORTER2_SETUP {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://jiarong/virsorter:2.2.3' :
        'jiarong/virsorter:2.2.3' }"

    input:
    val(meta)

    output:
    path "*"                                   , type: 'dir', maxDepth: 0, emit: db                                     
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ? task.ext.prefix+meta.id : meta.id
    def url      = 'https://osf.io/v46sc/download'
    """    
    virsorter \
        setup \
        -d ${prefix} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            virsorter --version | tail -n1 | sed 's/.*version //g')
    END_VERSIONS
    """
}