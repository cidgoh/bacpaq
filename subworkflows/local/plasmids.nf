// INCLUDE MODULES
include { MOBSUITE_RECON } from '../../modules/nf-core/mobsuite/recon/main'  
include { PLASMIDFINDER } from '../../modules/nf-core/plasmidfinder/main'

// Plasmid ID workflow
workflow PLASMIDS {
    take: genome
    main:
        
        // inititalize version channel
        ch_versions = Channel.empty()

        if (!params.skip_mobsuite) {
            MOBSUITE_RECON(genome)
            ch_versions = ch_versions.mix(MOBSUITE_RECON.out.versions)
            ch_mobsuite_fa = MOBSUITE_RECON.out.plasmids
        }
        
        if (!params.skip_plasmidfinder) {
            PLASMIDFINDER(genome)
            ch_versions = ch_versions.mix(PLASMIDFINDER.out.versions)
            ch_plasmidfinder_fa = PLASMIDFINDER.out.plasmid_seq
        }
        
    emit:
        versions = ch_versions
        plasmidfinder_fa = ch_plasmidfinder_fa
        mobsuite_fa = ch_mobsuite_fa
}