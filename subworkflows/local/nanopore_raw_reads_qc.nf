include { PORECHOP_ABI              } from '../../modules/nf-core/porechop/abi'
include { RASUSA as RASUSA_NANOPORE } from '../../modules/nf-core/rasusa'
include { NANOCOMP                  } from '../../modules/nf-core/nanocomp'
include { PYCOQC                    } from '../../modules/nf-core/pycoqc'

workflow NANOPORE_RAW_READS_QC {
    take:
    ch_merged_reads
    nanopore_summary_file

    main:
    ch_versions = Channel.empty()
    ch_sub_reads_qc = Channel.empty()
    ch_qc_reads = ch_merged_reads

    if (!params.skip_porechop) {
        PORECHOP_ABI(ch_merged_reads)
        ch_qc_reads = PORECHOP_ABI.out.reads
        ch_versions = ch_versions.mix(PORECHOP_ABI.out.versions)
    }
    if (!params.skip_subsampling) {

        // ch_coverages = Channel.fromList(params.depth_cut_off.split(',').collect { it.trim().toDouble() })

        if (!params.skip_porechop) {
            PORECHOP_ABI.out.reads
                .map { meta, fastq ->
                    tuple(meta, fastq, params.subsampling_genomesize)
                }
                .set { ch_sub_reads_qc }
        }
        else {
            ch_merged_reads
                .map { tuple(it[0], it[1], params.subsampling_genomesize) }
                .set { ch_sub_reads_qc }
        }
        RASUSA_NANOPORE(ch_sub_reads_qc, params.depth_cut_off)
        ch_qc_reads = RASUSA_NANOPORE.out.reads
        ch_versions = ch_versions.mix(RASUSA_NANOPORE.out.versions)
    }
    if (!params.skip_quality_report) {
        if (!params.skip_nanocomp) {
            ch_qc_reads
                .map { [it[1]] }
                .collect()
                .map { reads -> [[id: params.nanopore_summary_file_id], reads] }
                .set { ch_collected_reads }

            NANOCOMP(ch_collected_reads)
            ch_versions = ch_versions.mix(NANOCOMP.out.versions)
        }
        if (!params.skip_pycoqc) {
            PYCOQC(
                [[id: params.nanopore_summary_file_id], nanopore_summary_file]
            )
            ch_versions = ch_versions.mix(PYCOQC.out.versions)
        }
    }

    emit:
    versions     = ch_versions
    merged_reads = ch_qc_reads
}
