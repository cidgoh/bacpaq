include { ROARY } from '../../modules/nf-core/roary/main'
include { PIRATE } from '../../modules/nf-core/pirate/main'
include { PANAROO_RUN } from '../../modules/nf-core/panaroo/run/main'
include { PEPPAN } from '../../modules/local/PEPPAN/main'

workflow PANGENOME_ANALYSIS {
    ch_gff = [
        [ id:'GCA_010673125.1', single_end:false ],
        [
            "/home/$USER/scratch/hackathon/GCA_010673125.1.gff",
            "/home/$USER/scratch/hackathon/GCA_010798115.1.gff"
            
        ]
    ]

    if (params.reference_gff) {
        ch_reference_gff = params.reference_gff
    } else {
        ch_reference_gff = []
    }

    take:
    // should be ch_gff from prokka output
    // ch_gff

    main:

    ch_versions = Channel.empty()
    
    // skip setup, also add to nextflow.config
    if(!params.skip_roary) {
        ROARY(ch_gff)
        ch_versions = ch_versions.mix( ROARY.out.versions.first() )
    }

    if (!params.skip_pirate) {
        PIRATE(ch_gff)
        ch_versions = ch_versions.mix( PIRATE.out.versions.first() )
    }
    
    if (!params.skip_panaroo) {
        PANAROO_RUN(ch_gff)
        ch_versions = ch_versions.mix( PANAROO_RUN.out.versions.first() )
    }

    if (!params.skip_peppan) {
        PEPPAN(ch_gff, ch_reference_gff)
        ch_versions = ch_versions.mix( PEPPAN.out.versions.first() )
    }

    emit:
    versions = ch_versions
}