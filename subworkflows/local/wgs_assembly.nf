// import modules
include { SHOVILL } from '../../modules/nf-core/shovill'
include { RENAME_SHOVILL } from '../../modules/local/rename_shovill'
// include {  } from '../../modules/local/rename_shovill'

workflow WGS_ASSEMBLY {
    take: 
    reads

    main:
    reads.filter{ it[0].single_end == false } | SHOVILL
    
    // rename contig files
    RENAME_SHOVILL(SHOVILL.out.contigs, reads.map{ "fa" })

    emit:
    contigs = RENAME_SHOVILL.out
    corrections = SHOVILL.out.corrections
    shovill_log = SHOVILL.out.log
    raw_contigs = SHOVILL.out.raw_contigs
    gfa = SHOVILL.out.gfa
    versions = SHOVILL.out.versions

}