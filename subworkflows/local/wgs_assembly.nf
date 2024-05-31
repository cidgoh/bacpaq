// import modules
include { SHOVILL } from '../../modules/nf-core/shovill'
include { RENAME_CTG as RENAME_CTG_SHOVILL; RENAME_CTG as RENAME_CTG_DRAGONFLYE } from '../../modules/local/rename_ctg'
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
        versions = Channel.empty()
    //if (params.mode=="nanopore"){
        nanopore_reads | DRAGONFLYE
        RENAME_CTG_DRAGONFLYE(DRAGONFLYE.out.contigs, nanopore_reads.map{ "fa" })

        contigs = contigs.mix(RENAME_CTG_DRAGONFLYE.out)
        //corrections = corrections.mix(DRAGONFLYE.out.corrections)
        assembly_log = assembly_log.mix(DRAGONFLYE.out.log)
        raw_contigs = raw_contigs.mix(DRAGONFLYE.out.raw_contigs)
        gfa = gfa.mix(DRAGONFLYE.out.gfa)
        dragonflye_info = dragonflye_info.mix(DRAGONFLYE.out.txt)
        versions = versions.mix(DRAGONFLYE.out.versions)

//    }
//    else{
        illumina_reads.filter{ it[0].single_end == false } | SHOVILL

        // rename contig files
        RENAME_CTG_SHOVILL(SHOVILL.out.contigs, illumina_reads.map{ "fa" })
        contigs = RENAME_CTG_SHOVILL.out.mix(RENAME_CTG_SHOVILL.out)
        corrections = corrections.mix(SHOVILL.out.corrections)
        assembly_log = assembly_log.mix(SHOVILL.out.log)
        raw_contigs = raw_contigs.mix(SHOVILL.out.raw_contigs)
        gfa = gfa.mix(SHOVILL.out.gfa)
        versions = SHOVILL.out.versions

//    }

    emit:
    ch_contigs = contigs
    ch_corrections = corrections
    ch_asembly_log = assembly_log
    ch_raw_contigs = raw_contigs
    ch_gfa = gfa
    ch_dragonflye_info = dragonflye_info
    ch_versions = versions
}
