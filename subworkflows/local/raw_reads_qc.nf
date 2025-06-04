include { RASUSA                     } from '../../modules/nf-core/rasusa'
include { FASTQC as RAW_FASTQC       } from '../../modules/nf-core/fastqc'
include { FASTP                      } from '../../modules/nf-core/fastp'
include { TRIMMOMATIC                } from '../../modules/nf-core/trimmomatic'
include { TRIMGALORE                 } from '../../modules/nf-core/trimgalore'
include { FASTQC as TRIM_FASTQC      } from '../../modules/nf-core/fastqc'
include { MULTIQC as TRIM_MULTIQC    } from '../../modules/nf-core/multiqc'
include { MULTIQC as RAW_MULTIQC     } from '../../modules/nf-core/multiqc'
// include { CAT_FASTQ                  } from '../../modules/nf-core/cat/fastq/main'
// include { CAT_CAT                    } from '../../modules/nf-core/cat/cat/main'
include { CONFINDR                   } from '../../modules/local/confindr/confindr'
include { AGGREGATE_CONFINDR_RESULTS } from '../../modules/local/confindr/aggregate_confindr_results'
/*ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
*/
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RAW_READS_QC {
    take:
    ch_raw_reads


    main:
    ch_versions = Channel.empty()
    trimmomatic_report = Channel.empty()
    trimmomatic_log = Channel.empty()
    fastp_report = Channel.empty()
    raw_fastqc = Channel.empty()
    trim_fastqc = Channel.empty()


    //
    // Validate database paths
    //

    if (!params.skip_confindr) {
        if (params.confindr_db == null || !Utils.fileExists(params.confindr_db)) {
            log.error "Path to Confindr database is not valid"
            exit 1
            }
        }

    // use rasusa to randomly subsample sequencing reads
    if (!params.skip_subsampling) {
        // ch_genomesize = Channel.of(params.subsampling_genomesize)
        // ch_coverages = Channel.fromList(params.depth_cut_off.split(',').collect { it.trim().toDouble() })

        ch_raw_reads
            .map { tuple(it[0], it[1], params.subsampling_genomesize) }
            .set { ch_sub_reads_qc }

        RASUSA(ch_sub_reads_qc, params.depth_cut_off)
        ch_short_reads_fastqc = RASUSA.out.reads
        ch_short_reads_trim = RASUSA.out.reads
        ch_versions = ch_versions.mix(RASUSA.out.versions)
    }
    else {
        ch_short_reads_fastqc = ch_raw_reads
        ch_short_reads_trim = ch_raw_reads
    }

    // FASTQC check for raw reads
    RAW_FASTQC(
        ch_short_reads_fastqc
    )
    raw_fastqc = RAW_FASTQC.out.zip
    ch_versions = ch_versions.mix(RAW_FASTQC.out.versions)

    // Raw read trimming
    if (!params.skip_trimming) {
        // ch_adapter = Channel.of( file(params.adapter_fasta, checkIfExists: true) )
        // trim reads using fastp, trimommatic or trimgalore
        if (params.trim_tool == "fastp") {
            FASTP(
                ch_short_reads_trim,
                params.adapter_fasta != "null" ? file(params.adapter_fasta, checkIfExists: true) : ch_short_reads_trim.map { [] },
                params.discard_trimmed_pass,
                params.save_trimmed_fail,
                params.save_merged
            )
            ch_trimmed_reads = FASTP.out.reads
            fastp_report = FASTP.out.json
            ch_versions = ch_versions.mix(FASTP.out.versions)
        }
        else if (params.trim_tool == "trimommatic") {
            TRIMMOMATIC(ch_short_reads_trim)
            ch_trimmed_reads = TRIMMOMATIC.out.trimmed_reads
            trimmomatic_report = TRIMMOMATIC.out.summary
            trimmomatic_log = TRIMMOMATIC.out.out_log
            ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions)
        }
        else if (params.trim_tool == "trimgalore") {
            TRIMGALORE(ch_short_reads_trim)
            ch_trimmed_reads = TRIMGALORE.out.reads
            ch_versions = ch_versions.mix(TRIMGALORE.out.versions)
        }
        else {
            ch_trimmed_reads = ch_short_reads_trim
        }

        // FASTQC check for trimmed reads
        TRIM_FASTQC(
            ch_trimmed_reads
        )
        trim_fastqc = TRIM_FASTQC.out.zip
        ch_versions = ch_versions.mix(TRIM_FASTQC.out.versions)
    }
    else {
        ch_trimmed_reads = ch_short_reads_trim
    }

    if (!params.skip_confindr) {
        //Use confindr to detect contamination
        ch_confindr_results = CONFINDR(ch_trimmed_reads, params.confindr_db)
        confindr_csv = CONFINDR.out.csv
        confindr_log = CONFINDR.out.log
        confindr_report = CONFINDR.out.report
        confindr_rmlst = CONFINDR.out.rmlst
        ch_versions = ch_versions.mix(CONFINDR.out.versions)

        //Aggregate confindr results
        ch_confindr_results = Channel.empty()
        ch_confindr_results = ch_confindr_results.mix(CONFINDR.out.report.collect { it[1] })
        AGGREGATE_CONFINDR_RESULTS(ch_confindr_results)
        aggregate_confidr = AGGREGATE_CONFINDR_RESULTS.out.csv
    }

    emit:
    ch_trimmed_reads
    ch_versions
    raw_fastqc
    trim_fastqc
    fastp_report
    trimmomatic_report
    trimmomatic_log
    confindr_csv
    confindr_log
    confindr_report
    confindr_rmlst
    aggregate_confidr
}
