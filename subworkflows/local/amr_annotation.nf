// Importing the necessary modules for the workflow

include { ABRITAMR_RUN         } from '../../modules/nf-core/abritamr/run/main'
include { RGI_CARDANNOTATION   } from '../../modules/nf-core/rgi/cardannotation/main'
include { RGI_MAIN             } from '../../modules/nf-core/rgi/main/main'
include { AMRFINDERPLUS_UPDATE } from '../../modules/nf-core/amrfinderplus/update/main'
include { AMRFINDERPLUS_RUN    } from '../../modules/nf-core/amrfinderplus/run/main'
include { ABRICATE_RUN         } from '../../modules/nf-core/abricate/run/main'
include { ABRICATE_SUMMARY     } from '../../modules/nf-core/abricate/summary/main'
include { RESFINDER            } from '../../modules/local/resfinder/main'
include { UNTAR                } from '../../modules/nf-core/untar/main'

// This will be implemented in the future when the database
// information is ready for hamronization tool to parse

/*
include { HAMRONIZATION_SUMMARIZE } from '../../modules/nf-core/hamronization/summarize/main'
include { HAMRONIZATION_AMRFINDERPLUS } from '../../modules/nf-core/hamronization/amrfinderplus/main'
include { HAMRONIZATION_RGI } from '../../modules/nf-core/hamronization/rgi/main'
include { HAMRONIZATION_ABRICATE } from '../../modules/nf-core/hamronization/abricate/main'
*/

// Defining the workflow
workflow AMR_ANNOTATION {
    take:
    genome

    main:

    genome.view { meta, fasta -> "Debug: Genome meta: ${meta}, fasta: ${fasta}" }
    // Initializing empty channels for the output files
    ch_versions = Channel.empty()
    abricate_report = Channel.empty()
    abricate_summary_report = Channel.empty()
    amrfinderplus_report = Channel.empty()
    resfinder_report = Channel.empty()



    // Running abricate
    if (!params.skip_abricate) {
        ABRICATE_RUN(genome, [])
        abricate_report = ABRICATE_RUN.out.report
        ch_versions = ABRICATE_RUN.out.versions
        if (!params.skip_abricate_summary) {
            reports = abricate_report.collect()
            ABRICATE_SUMMARY(reports)
            abricate_summary_report = ABRICATE_SUMMARY.out.report
            ch_versions = ch_versions.mix(ABRICATE_SUMMARY.out.versions)
        }
    }
    rgi_tsv = Channel.empty()
    rgi_json = Channel.empty()
    rgi_db = Channel.empty()
    if (!params.skip_rgi) {
        card_db = [[id: 'card_db'], params.card_db]
        UNTAR(card_db)
        RGI_CARDANNOTATION(UNTAR.out.untar.map { it[1] })
        RGI_MAIN(genome, RGI_CARDANNOTATION.out.db, [])
        rgi_tsv = RGI_MAIN.out.tsv
        rgi_json = RGI_MAIN.out.json
        rgi_db = RGI_MAIN.out.tool_version
        ch_versions = ch_versions.mix(RGI_MAIN.out.versions)
    }

    if (!params.skip_amrfinderplus) {
        if (!params.amrfinderplus_db) {
            AMRFINDERPLUS_UPDATE()
            amrfinderplus_db = AMRFINDERPLUS_UPDATE.out.db
            ch_versions = ch_versions.mix(AMRFINDERPLUS_UPDATE.out.versions)
        }
        else {
            amrfinderplus_db = Channel.from(params.amrfinderplus_db)
        }
        AMRFINDERPLUS_RUN(genome, amrfinderplus_db)
        amrfinderplus_report = AMRFINDERPLUS_RUN.out.report
        ch_versions = ch_versions.mix(AMRFINDERPLUS_RUN.out.versions)
    }

    abritamr_matches = Channel.empty()
    abritamr_partials = Channel.empty()
    abritamr_virulence = Channel.empty()
    abritamr_amrfinder = Channel.empty()
    abritamr_txt = Channel.empty()
    if (!params.skip_abritamr) {
        ABRITAMR_RUN(genome)
        abritamr_matches = ABRITAMR_RUN.out.matches
        abritamr_partials = ABRITAMR_RUN.out.partials
        abritamr_virulence = ABRITAMR_RUN.out.virulence
        abritamr_amrfinder = ABRITAMR_RUN.out.out
        abritamr_txt = ABRITAMR_RUN.out.txt
        ch_versions = ch_versions.mix(ABRITAMR_RUN.out.versions)
    }

        if (!params.skip_resfinder){
            // Resfinder db validation
            if ( params.resfinder_db == null || !Utils.fileExists(params.resfinder_db)) {
                log.error "Path to ResFinder database was not provided or is not valid"
                exit 1
                }
            if ( params.pointfinder_db == null || !Utils.fileExists(params.pointfinder_db)) {
                log.error "Path to PointFinder database was not provided or is not valid"
                exit 1
                }

            RESFINDER(genome)
            resfinder_report = RESFINDER.out.resfinder_results_table
        }

    emit:
    abricate_report
    abricate_summary_report
    rgi_tsv
    rgi_json
    amrfinderplus_report
    abritamr_matches
    abritamr_partials
    abritamr_virulence
    abritamr_amrfinder
    abritamr_txt
    resfinder_report
    versions                = ch_versions
}
