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
    if (file ${barcode_dir}/*{fq,fastq}* | grep gzip) ; then 
        cat `find ${barcode_dir}/ -name "*fastq*" -or -name "*fq*"` > ${prefix}.merged.fastq.gz
    else 
        cat `find ${barcode_dir}/ -name "*fastq*" -or -name "*fq*"` | gzip > ${prefix}.merged.fastq.gz
    fi 
    """
}