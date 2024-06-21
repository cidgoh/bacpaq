// INCLUDING MODULES
include { RMLST } from '../../modules/local/RMLST'
include { RMLSTPARSE } from '../../modules/local/RMLSTPARSE'
//workflow

workflow RSMLST{
    take:
    assembly

    main:

    RMLST(assembly)
    RMLSTPARSE(RMLST.out.rmlstJSON)

    emit:
    json = RMLST.out.rmlstJSON
    summary = RMLSTPARSE.out.rmlstSUM
    allele = RMLSTPARSE.out.rmlstALLE
}
