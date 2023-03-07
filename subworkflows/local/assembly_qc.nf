// import modules
include { CHECKM_LINEAGEWF } from '../../modules/nf-core/checkm/lineagewf'
include { QUAST } from '../../modules/nf-core/quast'

workflow ASSEMBLY_QC {
    take: 
    assembly

    main:

    // RUN CHECKM
    assembly
        .map { file(it[1]).getExtension() }
        .set { assembly_ext }

    CHECKM_LINEAGEWF(
        assembly,
        assembly_ext,
        params.checkm_db
    )

    // RUN QUAST
    QUAST(
        assembly
    )

    emit:
    // CHECKM OUTPUTS
    checkm_output = CHECKM_LINEAGEWF.out.checkm_output
    marker_file = CHECKM_LINEAGEWF.out.marker_file
    checkm_tsv = CHECKM_LINEAGEWF.out.checkm_tsv
    versions = CHECKM_LINEAGEWF.out.versions    
    // QUAST OUTPUTS

}