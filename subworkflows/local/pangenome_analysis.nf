include { ROARY } from '../../modules/nf-core/roary/main'
include { PIRATE } from '../../modules/nf-core/pirate/main'

workflow PANGENOME_ANALYSIS {
    ch_gff = [
        [ id:'GCA_010673125.1', single_end:false ],
        [
            "/home/jyh25/scratch/hackathon/GCA_010673125.1.gff",
            "/home/jyh25/scratch/hackathon/GCA_010798115.1.gff"
        ]
    ]
    
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

    emit:
    // what happens if .out.results are null? best practices for these?
    // roary = ROARY.out.results
    // pirate = PIRATE.out.results
    versions = ch_versions
}