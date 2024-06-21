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
    ch_multiqc_files  = Channel.empty()
    ch_reads = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //


    ch_input
        .branch {
            illumina : it[0].mode == 'illumina'
            nanopore : it[0].mode == 'nanopore'
        }
        .set {ch_raw_reads}



    CAT_NANOPORE_FASTQ(ch_raw_reads.nanopore)
    ch_reads_merged = CAT_NANOPORE_FASTQ.out.reads
    //ch_versions = ch_versions.mix(CAT_NANOPORE_FASTQ.out.versions)
    ch_reads = ch_reads.mix(ch_raw_reads.illumina).mix(ch_reads_merged)

    //
    // SUBWORKFLOW: QC sub-workflow
    //
    if (!params.skip_QC){
        // SUBWORKFLOW: Run raw reads QC on nanopore & Illumina reads

        NANOPORE_RAW_READS_QC(
            ch_reads_merged,
            params.nanopore_summary_file
        )

        ch_versions = ch_versions.mix(NANOPORE_RAW_READS_QC.out.versions)

        RAW_READS_QC(ch_raw_reads.illumina)
        ch_tax_reads = ch_tax_reads.mix(RAW_READS_QC.out.short_reads)
        ch_tax_reads = ch_tax_reads.mix(NANOPORE_RAW_READS_QC.out.merged_reads)
        ch_assembly_reads = ch_tax_reads
        ch_versions = ch_versions.mix(RAW_READS_QC.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(RAW_READS_QC.out.raw_fastqc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(RAW_READS_QC.out.trim_fastqc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(RAW_READS_QC.out.fastp_report.collect{it[1]}.ifEmpty([]))
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
        ch_assembly_reads = TAXONOMY_QC.out.reads
        ch_versions = ch_versions.mix(TAXONOMY_QC.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(TAXONOMY_QC.out.kraken_report.collect{it[1]}.ifEmpty([]))
        //ch_multiqc_files = ch_multiqc_files.mix(TAXONOMY_QC.out.bracken_report.collect{it[1]}.ifEmpty([]))
    }
    /*else {
        ch_assembly_reads = ch_reads
    }*/

    //if (!params.skip_qc && !params.skip_taxonomy_qc) {
    ch_assembly_reads
        .branch {
        illumina : it[0].mode == 'illumina'
        nanopore : it[0].mode == 'nanopore'
        }
    .set {ch_assembly_reads}
    //}

    if (!params.skip_assembly) {
        // SUBWORKFLOW: Run WGS ASSEMBLY on reads
        WGS_ASSEMBLY(
            ch_assembly_reads.illumina,
            ch_assembly_reads.nanopore
        )
        ch_contigs = WGS_ASSEMBLY.out.ch_contigs
        ch_versions = ch_versions.mix(WGS_ASSEMBLY.out.versions)


        if (!params.skip_assembly_qc){
            // SUBWORKFLOW: Do ribosomal MLST on assembled contigs, using BIGSdb Restful API
            if (!params.skip_rmlst){
                RSMLST(
                    ch_contigs
                )
                //ch_versions = ch_versions.mix(RSMLST.out.versions)
            }
            // SUBWORKFLOW: RUN ASSEMBLY QC on assemblies
            ASSEMBLY_QC(
                ch_contigs
            )
            ch_versions = ch_versions.mix(ASSEMBLY_QC.out.versions)
        }
    }

    emit:
    contigs         = ch_contigs
    multiqc_files   = ch_multiqc_files
    versions        = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
