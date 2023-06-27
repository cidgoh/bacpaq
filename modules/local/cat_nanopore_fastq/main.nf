process CAT_NANOPORE_FASTQ {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(barcode_dir)

    output:
    tuple val(meta), path("*.merged.fastq.gz"), emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    """
    cat ${barcode_dir}/*.fastq.gz > ${prefix}.merged.fastq.gz
    """
}