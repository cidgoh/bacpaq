// INCLUDE MODULES
include { VIRSORTER2 } from '../../modules/local/virsorter2/main.nf'  

// Plasmid ID workflow
workflow PHAGE {
    take: genome
    main:
        
        // inititalize version channel
        ch_versions = Channel.empty()

        if (!params.skip_virsorter2) {
            VIRSORTER2(genome)
            ch_versions = ch_versions.mix(VIRSORTER2.out.versions)
            ch_virsorter2_fa = VIRSORTER2.out.fasta
        }
                
    emit:
        versions = ch_versions
        virsorter2_fa = ch_virsorter2_fa
}