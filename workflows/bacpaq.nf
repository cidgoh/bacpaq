
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



/*ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
//ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_bacpaq_pipeline'

include { SEQQC     } from '../subworkflows/local/seqqc'
include { ANNOTATION } from '../subworkflows/local/genome_annotation'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BACPAQ {

    take: samplesheet

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_illumina = samplesheet
        .filter { meta, fastq_1, fastq_2, fastq_dir, genome ->
            fastq_1 != null
        }
        .map {  meta, fastq_1, fastq_2, fastq_dir, genome ->
            if (!fastq_2) {
                return [ [ id:meta.id, single_end: true, mode:'illumina' ], [fastq_1 ]]
            } else {
                return [ [ id:meta.id, single_end: false, mode:'illumina' ], [fastq_1, fastq_2 ]]
            }
        }
    ch_onp = samplesheet
        .filter { meta, fastq_1, fastq_2, fastq_dir, genome ->
            fastq_dir != null
        }
        .map {  meta, fastq_1, fastq_2, fastq_dir, genome ->
                return [ [ id:meta.id, single_end: true, mode:'nanopore' ], [fastq_dir ]]
        }
    ch_genome = samplesheet
        .filter { meta, fastq_1, fastq_2, fastq_dir, genome ->
            genome != null
        }
        .map {  meta, fastq_1, fastq_2, fastq_dir, genome ->
                return [ [ id:meta.id ],[genome ]]
        }
    ch_reads = ch_illumina
        .concat(ch_onp)


    if ( !params.skip_seqqc ) {
        SEQQC(ch_reads)
        ch_genome = SEQQC.out.contigs
        ch_multiqc_files = ch_multiqc_files.mix(SEQQC.out.multiqc_files)
        ch_versions = ch_versions.mix(SEQQC.out.versions)
        }

    if ( !params.skip_annotation ) {
        if ( ch_genome ) {
            ANNOTATION(ch_genome)
            ch_versions = ch_versions.mix(ANNOTATION.out.versions)
        } else {
            log.error "${workflow.manifest.name}: No genomes available for annotation, exiting."
        }
    }


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'bacpaq_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //

    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.fromPath("${workflow.projectDir}/docs/images/nf-core-bacpaq_logo_light.png", checkIfExists: true)

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    /*if (!params.assembly_input) {

        if ( !params.skip_clipping && params.clip_tool == 'adapterremoval' ) {
            ch_multiqc_files = ch_multiqc_files.mix(ADAPTERREMOVAL_PE.out.settings.collect{it[1]}.ifEmpty([]))
            ch_multiqc_files = ch_multiqc_files.mix(ADAPTERREMOVAL_SE.out.settings.collect{it[1]}.ifEmpty([]))

        } else if ( !params.skip_clipping && params.clip_tool == 'fastp' )  {
            ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]}.ifEmpty([]))
        }

        if (!(params.keep_phix && params.skip_clipping && !(params.host_genome || params.host_fasta))) {
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.collect{it[1]}.ifEmpty([]))
        }

        if ( params.host_fasta || params.host_genome ) {
            ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_HOST_REMOVAL_ALIGN.out.log.collect{it[1]}.ifEmpty([]))
        }

        if(!params.keep_phix) {
            ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_PHIX_REMOVAL_ALIGN.out.log.collect{it[1]}.ifEmpty([]))
        }

    }

    ch_multiqc_files = ch_multiqc_files.mix(CENTRIFUGE_KREPORT.out.kreport.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2.out.report.collect{it[1]}.ifEmpty([]))

    if (!params.skip_quast){
        ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.report.collect().ifEmpty([]))

        if ( !params.skip_binning ) {
            ch_multiqc_files = ch_multiqc_files.mix(QUAST_BINS.out.dir.collect().ifEmpty([]))
        }
    }

    if ( !params.skip_binning || params.ancient_dna ) {
        ch_multiqc_files = ch_multiqc_files.mix(BINNING_PREPARATION.out.bowtie2_assembly_multiqc.collect().ifEmpty([]))
    }

    if (!params.skip_binning && !params.skip_prokka){
        ch_multiqc_files = ch_multiqc_files.mix(PROKKA.out.txt.collect{it[1]}.ifEmpty([]))
    }

    if (!params.skip_binning && !params.skip_binqc && params.binqc_tool == 'busco'){
        ch_multiqc_files = ch_multiqc_files.mix(BUSCO_QC.out.multiqc.collect().ifEmpty([]))
    }
    */

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
