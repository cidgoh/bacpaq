include { CAT_NANOPORE_FASTQ } from '../../modules/local/cat_nanopore_fastq'

workflow NANOPORE_RAW_READS_QC {
    take:
    ch_barcode_dirs

    main:
    ch_versions = Channel.empty()

    CAT_NANOPORE_FASTQ(ch_barcode_dirs)

    emit:
    barcode_reads = CAT_NANOPORE_FASTQ.out.reads
}