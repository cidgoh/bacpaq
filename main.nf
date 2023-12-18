#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/seqqc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/seqqc

    Website: https://nf-co.re/seqqc
    Slack  : https://nfcore.slack.com/channels/seqqc
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check if database filepaths exists
public static Boolean fileExists(filename) {
    return new File(filename).exists()
}

WorkflowMain.initialise(workflow, params, log)

if (!params.skip_kraken2 && !fileExists(params.kraken2_db)) {
    log.error "Path to Kraken2 database is not valid"
    exit 1
}
if ((params.classifier == 'centrifuge') && !fileExists(params.centrifuge_db)) {
    log.error "Path to Centrifuge database is not valid"
    exit 1
}
// if (!params.skip_checkm && !fileExists(params.checkm_db)) {
//     log.error "Path to CheckM database is not valid"
//     exit 1
// }
if (!params.skip_confindr && !fileExists(params.confindr_db)) {
    log.error "Path to Confindr database is not valid"
    exit 1
}
if (!params.skip_busco && !fileExists(params.busco_lineages_path)) {
    log.error "Path to BUSCO lineages database is not valid"
    exit 1
}

//include { validateParameters; paramsHelp; paramsSummaryLog; } from 'plugin/nf-validation'

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run my_pipeline --input input_file.csv")
   exit 0
}

// Validate input parameters
//validateParameters()

// Print summary of supplied parameters
//log.info paramsSummaryLog(workflow)

// Create a new channel of metadata from a sample sheet
// NB: `input` corresponds to `params.input` and associated sample sheet schema


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SEQQC } from './workflows/seqqc'
include { ANNOTATION } from './workflows/genome_annotation'
//
// WORKFLOW: Run main nf-core/seqqc analysis pipeline
//

workflow {

    main:  
    if (params.workflow == 'seqqc') {
        SEQQC ()
    }
    else if (params.workflow == 'annotation') {
        ANNOTATION ()
    }   
    else {
        log.error "Workflow not recognised"
        exit 1
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
/*workflow {
    NFCORE_SEQQC ()
}
*/
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
