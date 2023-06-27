include { PORECHOP_ABI } from '../../modules/nf-core/porechop/abi'
include { CAT_NANOPORE_FASTQ } from '../../modules/local/cat_nanopore_fastq'

workflow NANOPORE_RAW_READS_QC {
    take:
    ch_barcode_dirs

    main:
    ch_versions = Channel.empty()

    // Concatenate fastq files
    CAT_NANOPORE_FASTQ(ch_barcode_dirs)

    if (!params.skip_porechop) {
        PORECHOP_ABI(CAT_NANOPORE_FASTQ.out.reads)
    }

    emit:
    merged_reads = PORECHOP_ABI.out.reads
}