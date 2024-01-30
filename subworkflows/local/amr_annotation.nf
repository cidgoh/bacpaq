// Importing the necessary modules for the workflow

include { ABRITAMR_RUN } from '../../modules/nf-core/abritamr/run/main'
include { RGI_MAIN } from '../../modules/nf-core/rgi/main/main'
include { AMRFINDERPLUS_UPDATE } from '../../modules/nf-core/amrfinderplus/update/main'
include { AMRFINDERPLUS_RUN } from '../../modules/nf-core/amrfinderplus/run/main'
include { ABRICATE_RUN } from '../../modules/nf-core/abricate/run/main'
include { ABRICATE_SUMMARY } from '../../modules/nf-core/abricate/summary/main'

// This will be implemented in the future when the database
// information is ready for hamronization tool to parse

/*
include { HAMRONIZATION_SUMMARIZE } from '../../modules/nf-core/hamronization/summarize/main'
include { HAMRONIZATION_AMRFINDERPLUS } from '../../modules/nf-core/hamronization/amrfinderplus/main'
include { HAMRONIZATION_RGI } from '../../modules/nf-core/hamronization/rgi/main'
include { HAMRONIZATION_ABRICATE } from '../../modules/nf-core/hamronization/abricate/main'
*/

// Defining the workflow
workflow AMR_ANNOTATION{
    // Defining the input channel
    take:
        genome
    // Defining the main process
    main:

        // Initializing empty channels for the output files
        ch_versions = Channel.empty()
        abricate_report = Channel.empty()
        abricate_summary_report = Channel.empty()
        rgi_report = Channel.empty()
        amrfinderplus_report = Channel.empty()
        abritamr_report = Channel.empty()


        // Running abricate
        if (!params.skip_abricate){
            abricate_dbs = Channel.empty()
            ABRICATE_RUN(genome)
            abricate_report = ABRICATE_RUN.out.report
            ch_versions = ABRICATE_RUN.out.versions
            if (!params.skip_abricate_summary){
                reports = abricate_report.collect()
                ABRICATE_SUMMARY(reports)
                abricate_summary_report = ABRICATE_SUMMARY.out.report
                ch_versions = ch_versions.mix(ABRICATE_SUMMARY.out.versions)
            }
            /*
            if (!params.skip_hamronization ){
                HAMRONIZATION_ABRICATE(abricate_summary_report, "tsv", ABRICATE_RUN.out.versions, "test")
            }*/

        }

        if (!params.skip_rgi){
            RGI_MAIN(genome)
            rgi_report = RGI_MAIN.out.tsv
            ch_versions = ch_versions.mix(RGI_MAIN.out.versions)

            if (!params.skip_hamronization ){
                HAMRONIZATION_RGI(rgi_report, "tsv", RGI_MAIN.out.versions, "test")
            }
        }

        if (!params.skip_amrfinderplus){
            if (!params.amrfinderplus_db){
                AMRFINDERPLUS_UPDATE()
                amrfinderplus_db = AMRFINDERPLUS_UPDATE.out.db
                ch_versions = ch_versions.mix(AMRFINDERPLUS_UPDATE.out.versions)
            }
            else{
                amrfinderplus_db = Channel.from(params.amrfinderplus_db)

            }
            AMRFINDERPLUS_RUN(genome, amrfinderplus_db)
            amrfinderplus_report = AMRFINDERPLUS_RUN.out.report
            ch_versions = ch_versions.mix(AMRFINDERPLUS_RUN.out.versions)


            if (!params.skip_hamronization ){
                HAMRONIZATION_AMRFINDERPLUS(amrfinderplus_report, "tsv", "test", "test")
            }

        }

        if (!params.skip_abritamr){
            ABRITAMR_RUN(genome)
            abritamr_report = ABRITAMR_RUN.out.txt
            ch_versions = ch_versions.mix(ABRITAMR_RUN.out.versions)
        }




    // Defining the output channels
    emit:
        abricate_report
        abricate_summary_report
        rgi_report
        amrfinderplus_report
        abritamr_report
        versions = ch_versions
}
