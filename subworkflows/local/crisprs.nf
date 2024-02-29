include { CCTYPER } from '../../modules/local/cctyper/main'

workflow CRISPRS {
    take: genome
    main:
        // initialize versions channel
        ch_versions = Channel.empty()

        // RUN CCTYPER
        if (!params.skip_cctyper) {
            CCTYPER(genome)
            ch_versions = ch_versions.mix(CCTYPER.out.versions)
        }

    emit:
        gff = CCTYPER.out.gff
        versions = ch_versions
}
