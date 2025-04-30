// import modules
include { SHOVILL    } from '../../modules/nf-core/shovill'
include { RENAME_CTG } from '../../modules/local/rename_ctg'
// include {  } from '../../modules/local/rename_shovill'
include { DRAGONFLYE } from '../../modules/nf-core/dragonflye'

workflow WGS_ASSEMBLY {
    take:
    illumina_reads
    nanopore_reads

    main:

    contigs = Channel.empty()
    corrections = Channel.empty()
    assembly_log = Channel.empty()
    raw_contigs = Channel.empty()
    gfa = Channel.empty()
    dragonflye_info = Channel.empty()
    ch_versions = Channel.empty()
    renamed_contigs = Channel.empty()

    //if (params.mode=="nanopore"){

    DRAGONFLYE(
        nanopore_reads.map { meta, file ->
            [meta, [], file]
        }
    )
    illumina_reads.filter { it[0].single_end == false } | SHOVILL
    contigs = contigs.mix(DRAGONFLYE.out.contigs)
    contigs = contigs.mix(SHOVILL.out.contigs)
    RENAME_CTG(
        contigs.map { tuple(it[0], it[1]) },
        contigs.map { it[1].getName().split('\\.')[-1] },
    )
    //RENAME_CTG(contigs, contigs.map { it[1].getExtension() })
    renamed_contigs = renamed_contigs.mix(RENAME_CTG.out)

    //corrections = corrections.mix(DRAGONFLYE.out.corrections)
    assembly_log = assembly_log.mix(DRAGONFLYE.out.log)
    raw_contigs = raw_contigs.mix(DRAGONFLYE.out.raw_contigs)
    gfa = gfa.mix(DRAGONFLYE.out.gfa)
    dragonflye_info = dragonflye_info.mix(DRAGONFLYE.out.txt)
    ch_versions = ch_versions.mix(DRAGONFLYE.out.versions)


    corrections = corrections.mix(SHOVILL.out.corrections)
    assembly_log = assembly_log.mix(SHOVILL.out.log)
    raw_contigs = raw_contigs.mix(SHOVILL.out.raw_contigs)
    gfa = gfa.mix(SHOVILL.out.gfa)
    ch_versions = ch_versions.mix(SHOVILL.out.versions)

    emit:
    ch_contigs         = renamed_contigs
    ch_corrections     = corrections
    ch_asembly_log     = assembly_log
    ch_raw_contigs     = raw_contigs
    ch_gfa             = gfa
    ch_dragonflye_info = dragonflye_info
    versions           = ch_versions
}
