include { PORECHOP_ABI } from '../../modules/nf-core/porechop/abi'
include { CAT_NANOPORE_FASTQ } from '../../modules/local/cat_nanopore_fastq'
include { RASUSA as RASUSA_NANOPORE } from '../../modules/nf-core/rasusa'
include { NANOCOMP } from '../../modules/nf-core/nanocomp'

workflow NANOPORE_RAW_READS_QC {
    take:
    ch_barcode_dirs

    main:
    ch_versions = Channel.empty()

    // Concatenate fastq files
    CAT_NANOPORE_FASTQ(ch_barcode_dirs)

    ch_merged_reads = CAT_NANOPORE_FASTQ.out.reads

    if (!params.skip_porechop) {
        PORECHOP_ABI(CAT_NANOPORE_FASTQ.out.reads)
        ch_merged_reads = PORECHOP_ABI.out.reads
    }
    if (!params.skip_subsampling) {
        ch_genomesize = Channel.of(params.genomesize)

        ch_merged_reads
                    .combine(ch_genomesize)
                    .set { ch_sub_reads_qc }

        RASUSA_NANOPORE(ch_sub_reads_qc, params.depth_cut_off)

        ch_merged_reads = RASUSA_NANOPORE.out.reads
    }
    if (!params.skip_quality_report) {
        NANOCOMP(ch_merged_reads)
    }

    emit:
    merged_reads = ch_merged_reads
}