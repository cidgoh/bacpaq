// import modules
include { SHOVILL } from '../../modules/nf-core/shovill'

workflow WGS_ASSEMBLY {
    take: 
    reads

    main:
    reads.filter{ it[0].single_end == false } | SHOVILL

    emit:
    contigs = SHOVILL.out.contigs
    corrections = SHOVILL.out.corrections
    shovill_log = SHOVILL.out.log
    raw_contigs = SHOVILL.out.raw_contigs
    gfa = SHOVILL.out.gfa
    versions = SHOVILL.out.versions

}