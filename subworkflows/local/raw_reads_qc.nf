include { FASTP                        } from '../../modules/nf-core/fastp'
include { FASTQC                    } from '../../modules/nf-core/fastqc'
include { MULTIQC                           } from '../../modules/nf-core/multiqc'
include { CAT_FASTQ                   } from '../../modules/nf-core/cat/fastq/main'
include { CAT_CAT } from '../../modules/nf-core/cat/cat/main' 


workflow RAW_READS_QC {
    take:
    ch_raw_reads_qc

    main:
    ch_versions = Channel.empty()
    if (params.trimming_tool=="fastp") {

        ch_adaptor=Channel.from(params.adapter_fasta)
        FASTP (
            ch_raw_reads_qc, ch_adaptor, params.save_trimmed_fail, params.save_merged
        )
        // ch_short_reads = FASTP.out.reads
        // ch_versions = ch_versions.mix(FASTP.out.versions.first())
    }
    else {

    }

    // INPUT_CHECK (
    //     ch_input
    // )
    // .reads
    // .map {
    //     meta, fastq ->
    //         new_id = meta.id - ~/_T\d+/
    //         [ meta + [id: new_id], fastq ]
    // }
    // .groupTuple()
    // .branch {
    //     meta, fastq ->
    //         single  : fastq.size() == 1
    //             return [ meta, fastq.flatten() ]
    //         multiple: fastq.size() > 1
    //             return [ meta, fastq.flatten() ]
    // }
    // .set { ch_fastq }
    // ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // ch_fastq.view()
    // //
    // // MODULE: Concatenate FastQ files from same sample if required
    // //
    // CAT_FASTQ (
    //     ch_fastq.multiple
    // )
    // .reads
    // .mix(ch_fastq.single)
    // .set { ch_cat_fastq }
    // ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))
    

}


