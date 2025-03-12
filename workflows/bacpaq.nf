
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
def convertEmptyToNull = { value ->
    if (value.toString().trim().isEmpty() || (value instanceof List && value.isEmpty())) {
        return null
    }
    return value
}

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
include { VARIANT_CALLING } from '../subworkflows/local/variantcalling'
include { VARIANT_VIS } from '../subworkflows/local/variantvis'
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
        .map { meta, fastq_1, fastq_2, fastq_dir, genome ->
            [meta, convertEmptyToNull(fastq_1), fastq_2, fastq_dir, genome]
        }
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
    ch_illumina.ifEmpty {
        log.warn "No Illumina samples found in the samplesheet, skipping Illumina-specifc processes."
    }
    ch_onp = samplesheet
        .map { meta, fastq_1, fastq_2, fastq_dir, genome ->
            [meta, fastq_1, fastq_2, convertEmptyToNull(fastq_dir), genome]
        }
        .filter { meta, fastq_1, fastq_2, fastq_dir, genome ->
            fastq_dir != null
        }
        .map {  meta, fastq_1, fastq_2, fastq_dir, genome ->
                return [ [ id:meta.id, single_end: true, mode:'nanopore' ], [fastq_dir ]]
        }
    ch_onp.ifEmpty {
        log.warn "No Nanopore samples found in the samplesheet, skipping Nanopore-specifc processes."
    }
    ch_genome = samplesheet
        .map { meta, fastq_1, fastq_2, fastq_dir, genome ->
            [meta, fastq_1, fastq_2, fastq_dir, convertEmptyToNull(genome)]
        }
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
            ch_multiqc_files = ch_multiqc_files.mix(ANNOTATION.out.multiqc_files)
            ch_versions = ch_versions.mix(ANNOTATION.out.versions)
        } else {
            log.error "${workflow.manifest.name}: No genomes available for annotation, exiting."
        }
    }

    //
    // VARIANT DETECTION WORKFLOW
    //
    if (!params.skip_variant_detection) {
        // collect all inputs into a single channel
        ch_variant_input = ch_genome
            .concat(ch_onp)
            .concat(ch_illumina)
            
        // VARIANT CALLING SUBWORKFLOW
        VARIANT_CALLING(ch_variant_input)

        // VARIANT VIZ SUBWORKFLOW
        if (!params.skip_variant_viz) {

            ch_vcf = VARIANT_CALLING.out.vcf_snippy
                .concat(VARIANT_CALLING.out.vcf_medaka)
                .concat(VARIANT_CALLING.out.vcf_nucmer)
            ch_bam = VARIANT_CALLING.out.bam_snippy
                .concat(VARIANT_CALLING.out.bam_medaka)
                .concat(VARIANT_CALLING.out.bam_nucmer_sorted)
            ch_bai = VARIANT_CALLING.out.bai_snippy
                .concat(VARIANT_CALLING.out.bai_medaka)
                .concat(VARIANT_CALLING.out.bai_nucmer)
            ch_aln_fa = VARIANT_CALLING.out.core_aln

            VARIANT_VIS(
                ch_vcf,
                ch_bam,
                ch_bai,
                ch_aln_fa
            )
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
        Channel.fromPath("${workflow.projectDir}/docs/images/LogoCIDGOH2.png", checkIfExists: true)

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
