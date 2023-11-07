// INCLUDING MODULES
// Importing the necessary modules for the workflow
 
include { ABRITAMR_RUN } from '../../modules/nf-core/abritamr/run/main'                                                                     
include { RGI_MAIN } from '../../modules/nf-core/rgi/main/main'
include { AMRFINDERPLUS_UPDATE } from '../../modules/nf-core/amrfinderplus/update/main'
include { AMRFINDERPLUS_RUN } from '../../modules/nf-core/amrfinderplus/run/main'
include { ARIBA_GETREF } from '../../modules/nf-core/ariba/getref/main'
include { ARIBA_RUN } from '../../modules/nf-core/ariba/run/main'
include { ABRICATE_RUN } from '../../modules/nf-core/abricate/run/main'
include { ABRICATE_SUMMARY } from '../../modules/nf-core/abricate/summary/main'

// Defining the workflow
workflow AMR_ANNOTATION{
    // Defining the input channel
    take: 
        ch_genome

    // Defining the main process
    main:
        // Initializing empty channels for the output files
        ch_versions = Channel.empty()
        abricate_report = Channel.empty()

        // Running abricate
        if (!params.skip_abricate){
            ABRICATE_RUN(ch_genome)
            abricate_report = ABRICATE_RUN.out.report
            ch_versions = ABRICATE_RUN.out.versions
            if (!params.skip_abricate_summary){
                reports = abricate_report.collect()
                ABRICATE_SUMMARY(reports)
                summarized_report = ABRICATE_SUMMARY.out.report
                ch_versions = ch_versions.mix(ABRICATE_SUMMARY.out.versions)
            }
        }
        
        if (!params.skip_rgi){
            RGI_MAIN(ch_genome)
            rgi_report = RGI_MAIN.out.tsv
            ch_versions = ch_versions.mix(RGI_MAIN.out.versions)
        }


        
    // Defining the output channels
    emit:
        abricate_report
        ch_versions
}