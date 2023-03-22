include { FASTP                    } from '../../modules/nf-core/fastp'
include { TRIMMOMATIC              } from '../../modules/nf-core/trimmomatic'
include { TRIMGALORE               } from '../../modules/nf-core/trimgalore'
include { FASTQC as RAW_FASTQC     } from '../../modules/nf-core/fastqc'
include { FASTQC as TRIM_FASTQC    } from '../../modules/nf-core/fastqc'
include { MULTIQC as TRIM_MULTIQC  } from '../../modules/nf-core/multiqc'
include { MULTIQC as RAW_MULTIQC  } from '../../modules/nf-core/multiqc'
include { RASUSA                   } from '../../modules/nf-core/rasusa'
include { CONFINDR                 } from '../../modules/local/confindr'
include { AGGREGATE_CONFINDR_RESULTS                } from '../../modules/local/aggregate_confindr_results'
include { CAT_FASTQ                } from '../../modules/nf-core/cat/fastq/main'
include { CAT_CAT                  } from '../../modules/nf-core/cat/cat/main'

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)


workflow RAW_READS_QC {
    take:
    ch_raw_reads_qc

    main:
    ch_versions = Channel.empty()

    // use rasusa to randomly subsample sequencing reads
    if (!params.skip_subsampling) {

        ch_genomesize= Channel.of(params.genomesize)

        ch_raw_reads_qc
                    .combine(ch_genomesize)
                    .set { ch_sub_reads_qc }

        RASUSA (ch_sub_reads_qc, params.depth_cut_off)

        ch_raw_reads_qc = RASUSA.out.reads
    }

    // FASTQC check for raw reads
    RAW_FASTQC (
        ch_raw_reads_qc
    )
    ch_versions = ch_versions.mix(RAW_FASTQC.out.versions.first())

    // trim reads using fastp, trimommatic or trimgalore
    if (params.trim_tool=="fastp") {

        adaptor=file params.adapter_fasta
        FASTP (
            ch_raw_reads_qc, adaptor, params.save_trimmed_fail, params.save_merged
        )
        ch_short_reads = FASTP.out.reads
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
    }
    else if (params.trim_tool=="trimommatic") {

        TRIMMOMATIC (ch_raw_reads_qc)
        ch_short_reads = TRIMMOMATIC.out.trimmed_reads
        ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions.first())
    }
    else if (params.trim_tool=="trimgalore") {

        TRIMGALORE (ch_raw_reads_qc)
        ch_short_reads = TRIMGALORE.out.reads
        ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())
    }

    // FASTQC check for trimmed reads
    TRIM_FASTQC (
        ch_short_reads
    )
    ch_versions = ch_versions.mix(TRIM_FASTQC.out.versions.first())

    if(!params.skip_confindr){
        //Use confindr to detect contamination
        ch_confindr_results = CONFINDR(ch_raw_reads_qc, params.confindr_db)
        ch_versions = ch_versions.mix(CONFINDR.out.versions.first())

        //Aggregate confindr results
        ch_confindr_results = Channel.empty()
        ch_confindr_results = ch_confindr_results.mix(CONFINDR.out.report.collect({it[1]}))
        AGGREGATE_CONFINDR_RESULTS(ch_confindr_results)
    }
    

    // MultiQC report for raw reads
    ch_raw_multiqc_files = Channel.empty()
    ch_raw_multiqc_files = ch_raw_multiqc_files.mix(RAW_FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    RAW_MULTIQC (
        ch_raw_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    // MultiQC report for trimmed reads
    ch_trim_multiqc_files = Channel.empty()
    ch_trim_multiqc_files = ch_trim_multiqc_files.mix(TRIM_FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    TRIM_MULTIQC (
        ch_trim_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    emit:
    short_reads=ch_short_reads


}


