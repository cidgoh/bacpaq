

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


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
    multiqc_report = ''


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
    ch_reads.view()

    if ( !params.skip_seqqc ) {
        SEQQC(ch_reads)
        ch_genome = SEQQC.out.ch_contigs
    }

    if ( !params.skip_annotation ) {
        if ( ch_genome ) {
            ANNOTATION(ch_genome)
        } else {
            log.error "${workflow.manifest.name}: No genomes available for annotation, exiting."
        }
    }
    // if (params.workflow == 'seqqc') {
    // }
    // else if (params.workflow == 'annotation') {

    //     samplesheet.view()
    // }
    // else {
    //     log.error "Workflow not recognised"
    //     exit 1
    // }

    // SEQQC ()
    // ANNOTATION ()


    emit:
    // Output files
    multiqc_report

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/*workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
