// INCLUDING MODULES
include { BAKTA_BAKTA } from '../../modules/nf-core/bakta/bakta/main' 
include { BAKTA_BAKTADBDOWNLOAD } from '../../modules/nf-core/bakta/baktadbdownload/main'
include { PROKKA } from '../../modules/nf-core/prokka/main'                                                                          
//workflow


workflow GENE_ANNOTATION{
    take: 
        ch_genome

    main:
        PROKKA(ch_genome, [], [])
        gff = PROKKA.out.gff
        //ch_versions   = ch_versions.mix(PROKKA.out.versions)

    emit:
        gff
        //ch_versions
    }
