//
// Importing the modules required for the sub-workflow
//


include { KRAKEN2_KRAKEN2 } from '../../modules/nf-core/kraken2/kraken2/main'
include { KRAKENTOOLS_COMBINEKREPORTS } from '../../modules/nf-core/krakentools/combinekreports/main'
include { KRAKENTOOLS_KREPORT2KRONA } from '../../modules/nf-core/krakentools/kreport2krona/main'
include { BRACKEN_BRACKEN } from '../../modules/nf-core/bracken/bracken/main'
include { BRACKEN_COMBINEBRACKENOUTPUTS } from '../../modules/nf-core/bracken/combinebrackenoutputs/main'
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main'
include { CENTRIFUGE_CENTRIFUGE } from '../../modules/nf-core/centrifuge/centrifuge/main'
include { CENTRIFUGE_KREPORT } from '../../modules/nf-core/centrifuge/kreport/main'


workflow TAXONOMY_QC {
    take:
    ch_reads_taxonomy
    classified_reads
    unclassified_reads
    

    main:
    ch_db=Channel.from(params.kraken2_db)
    ch_versions = Channel.empty()
    
    if (!params.skip_kraken2) {
        KRAKEN2_KRAKEN2(
            ch_reads_taxonomy,
            ch_db,
            classified_reads,
            unclassified_reads
        )

    }

    if (!params.skip_combinekreports) {
        KRAKENTOOLS_COMBINEKREPORTS(
            KRAKEN2_KRAKEN2.out.report.collect())
            
    }  

    if (!params.skip_kreport2krona) {
        KRAKENTOOLS_KREPORT2KRONA(
            KRAKENTOOLS_COMBINEKREPORTS.out.txt)
            
    } 
    
    /*if (!params.skip_bracken) {
        BRACKEN_BRACKEN(
            KRAKEN2_KRAKEN2.out.report)
            
    } */


    //if (!params.skip_dehosting){

    //}
    
    emit:
    classified_reads                           // channel: [ val(meta), [ reads ] ]
    unclassified_reads                         // channel: [ val(meta), [ reads ] ]
    //versions = TAXONOMY_QC.out.versions.first() // channel: [ versions.yml ]
}
