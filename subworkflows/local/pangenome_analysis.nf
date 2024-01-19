include { ROARY } from '../../modules/nf-core/roary/main'
include { PIRATE } from '../../modules/nf-core/pirate/main'
include { PANAROO_RUN } from '../../modules/nf-core/panaroo/run/main'
include { PEPPAN } from '../../modules/local/PEPPAN/main'

workflow PANGENOME_ANALYSIS {


    take:
    // should be ch_gff from prokka output
    //ch_gff

    main:

    ch_versions = Channel.empty()

    ch_gff = [
        [ id:'pan_genome', single_end:false ],
        [
            file("/scratch/mzanwar/nf-seqqc/PROKKA_GCA_010673125.1/PROKKA_GCA_010673125.1.gff"),
            file("/scratch/mzanwar/nf-seqqc/PROKKA_GCA_010798115.1/PROKKA_GCA_010798115.1.gff")

        ]
    ]

    if (params.reference_gff) {
        ch_reference_gff = params.reference_gff
    } else {
        ch_reference_gff = []
    }

    // skip setup, also add to nextflow.config
    if(!params.skip_roary) {
        ROARY(ch_gff)
        roary_results = ROARY.out.results
        roary_alignment = ROARY.out.aln
        ch_versions = ch_versions.mix( ROARY.out.versions )
    }

    if (!params.skip_pirate) {
        PIRATE(ch_gff)
        pirate_results = PIRATE.out.results
        pirate_alignment = PIRATE.out.aln
        ch_versions = ch_versions.mix( PIRATE.out.versions )
    }

    if (!params.skip_panaroo) {
        PANAROO_RUN(ch_gff)
        panaroo_results = PANAROO_RUN.out.results
        panaroo_alignment = PANAROO_RUN.out.aln
        ch_versions = ch_versions.mix( PANAROO_RUN.out.versions )
    }

    if (!params.skip_peppan) {
        PEPPAN(ch_gff, ch_reference_gff)
        peppan_gff = PEPPAN.out.gff_peppan
        peppan_fna = PEPPAN.out.fna_peppan
        ch_versions = ch_versions.mix( PEPPAN.out.versions )
    }

    emit:
    versions = ch_versions
    roary_results
    roary_alignment
    pirate_results
    pirate_alignment
    panaroo_results
    panaroo_alignment
    peppan_gff
    peppan_fna
}
