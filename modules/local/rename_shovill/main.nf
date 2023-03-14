process RENAME_SHOVILL {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(assembly)
    val file_ext

    output:
    tuple val(meta), path("*")

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    prefix       = task.ext.prefix ?: "${meta.id}"
    """
    mv \
        $args \
        ${assembly} \
        ${prefix}.${file_ext}
    """
}
