/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check if database filepaths exists

//def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
//WorkflowSeqqc.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
// def checkPathParamList = [  params.multiqc_config, params.fasta ]
// for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }


// Check if input database paths are valid
if (!params.skip_kraken2 && !Utils.fileExists(params.kraken2_db)) {
    log.error "Path to Kraken2 database is not valid"
    exit 1
}
if ((params.classifier == 'centrifuge') && !Utils.fileExists(params.centrifuge_db)) {
    log.error "Path to Centrifuge database is not valid"
    exit 1
}
if (!params.skip_checkm && !Utils.fileExists(params.checkm_db)) {
    log.error "Path to CheckM database is not valid"
    exit 1
}
if (!params.skip_confindr && !Utils.fileExists(params.confindr_db)) {
    log.error "Path to Confindr database is not valid"
    exit 1
}
if (!params.skip_busco && !Utils.fileExists(params.busco_lineages_path)) {
    log.error "Path to BUSCO lineages database is not valid"
    exit 1
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
//include { INPUT_CHECK           } from './input_check'
include { WGS_ASSEMBLY          } from './wgs_assembly'
include { ASSEMBLY_QC           } from './assembly_qc'
include { RSMLST                } from './rmlst'
include { TAXONOMY_QC           } from './taxonomy_qc'
include { RAW_READS_QC          } from './raw_reads_qc'
include { CAT_NANOPORE_FASTQ    } from '../../modules/local/cat_nanopore_fastq/main'
include { NANOPORE_RAW_READS_QC } from './nanopore_raw_reads_qc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
//include { FASTQC                      } from '../../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../../modules/nf-core/multiqc/main'
//include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../../modules/nf-core/custom/dumpsoftwareversions/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SEQQC {

    take:
        ch_input

    main:
    ch_versions = Channel.empty()
    ch_tax_reads = Channel.empty()
    ch_assembly_reads = Channel.empty()
    ch_contigs = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    //ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    ch_input
        .branch {
            illumina : it[0].mode == 'illumina'
            nanopore : it[0].mode == 'nanopore'
        }
        .set {ch_reads}



    CAT_NANOPORE_FASTQ(ch_reads.nanopore)
    ch_reads_merged = CAT_NANOPORE_FASTQ.out.reads


    //
    // SUBWORKFLOW: QC sub-workflow
    //
    if (!params.skip_QC){
        // SUBWORKFLOW: Run raw reads QC on nanopore & Illumina reads

        NANOPORE_RAW_READS_QC(
            ch_reads_merged,
            params.nanopore_summary_file
        )

        RAW_READS_QC(ch_reads.illumina)
        ch_tax_reads = RAW_READS_QC.out.short_reads.mix(NANOPORE_RAW_READS_QC.out.merged_reads)
        ch_assembly_reads = ch_tax_reads
    }
    else{
        ch_tax_reads = ch_reads
        ch_assembly_reads = ch_reads
    }

    if(!params.skip_taxonomy_qc) {
        TAXONOMY_QC (
            ch_tax_reads,
            params.host_genome
            )
        ch_assembly_reads = TAXONOMY_QC.out
    } else {
        ch_assembly_reads = ch_reads
    }

    if (!params.skip_qc && !params.skip_taxonomy_qc) {
        ch_assembly_reads
            .branch {
            illumina : it[0].mode == 'illumina'
            nanopore : it[0].mode == 'nanopore'
            }
        .set {ch_assembly_reads}
    }

    if (!params.skip_assembly) {
        // SUBWORKFLOW: Run WGS ASSEMBLY on reads
        WGS_ASSEMBLY(
            ch_assembly_reads.illumina,
            ch_assembly_reads.nanopore
        )
        ch_contigs = WGS_ASSEMBLY.out.ch_contigs

        if (!params.skip_assembly_qc){
            // SUBWORKFLOW: Do ribosomal MLST on assembled contigs, using BIGSdb Restful API
            if (!params.skip_rmlst){
                RSMLST(
                    ch_contigs
                )
            }
            // SUBWORKFLOW: RUN ASSEMBLY QC on assemblies
            ASSEMBLY_QC(
                ch_contigs
            )
        }
    }
    //ch_versions = ch_versions.mix(TAXONOMY_QC.out.versions.first())

    //CUSTOM_DUMPSOFTWAREVERSIONS (
    //    ch_versions.unique().collectFile(name: 'collated_versions.yml')
    //)


    //
    // MODULE: MultiQC
    //
    /*workflow_summary    = WorkflowSeqqc.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowSeqqc.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    //ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report
    */

    emit:
    ch_contigs
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/* workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
} */

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
