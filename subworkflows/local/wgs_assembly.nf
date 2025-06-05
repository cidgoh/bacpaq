// import modules
include { SHOVILL    } from '../../modules/nf-core/shovill'
include { RENAME_CTG } from '../../modules/local/rename_ctg'
include { FLYE } from '../../modules/nf-core/flye'
include { GUNZIP; GUNZIP as GUNZIP_MEDAKA } from '../../modules/nf-core/gunzip'
include { MEDAKA } from '../../modules/nf-core/medaka/'

workflow WGS_ASSEMBLY {
    take:
    illumina_reads
    nanopore_reads

    main:

    corrections = Channel.empty()
    assembly_log = Channel.empty()
    raw_contigs = Channel.empty()
    gfa = Channel.empty()
    flye_info = Channel.empty()
    ch_versions = Channel.empty()
    renamed_contigs = Channel.empty() // polished contigs for both modes
    illumina_contigs = Channel.empty() // polished contigs for illumina
    nanopore_contigs = Channel.empty() // polished contigs for nanopore

    // RUN FLYE ASSEMBLY ON NANOPORE READS
    nanopore_qual = Channel.of("--nano-raw", "--nano-corr", "--nano-hq")
    nanopore_qual.filter(val -> val =~ params.nanopore_qual) // determine nanopore read quality [raw, corr, hq]
         .set { ch_nanopore_qual }
    FLYE(nanopore_reads, ch_nanopore_qual.first())
    GUNZIP(FLYE.out.fasta)
    // NANOPORE POLISHING WITH MEDAKA
    if (!params.skip_medaka) {
        medaka_repo = "https://github.com/nanoporetech/medaka/raw/refs/heads/master/medaka/data/"
        ch_medaka_model = Channel.from(params.medaka_model)
            .map { model ->
                medaka_repo+model+"_model.tar.gz"
            }
        ch_medaka = nanopore_reads.join(GUNZIP.out.gunzip)
        MEDAKA(
            ch_medaka, 
            ch_medaka_model.first()
        )
        GUNZIP_MEDAKA(
            MEDAKA.out.assembly
        )
        nanopore_contigs = nanopore_contigs.mix(GUNZIP_MEDAKA.out.gunzip)
    } else {
        nanopore_contigs = nanopore_contigs.mix(FLYE.out.fasta)
    }

    // RUN ILLUMINA READ ASSEMBLY
    illumina_reads.filter { it[0].single_end == false } | SHOVILL
    illumina_contigs = illumina_contigs.mix(SHOVILL.out.contigs) // polished contigs

    // rename contig file names for illumina assemblies ONLY
    RENAME_CTG(
        illumina_contigs.map { tuple(it[0], it[1]) },
        illumina_contigs.map { it[1].getName().split('\\.')[-1] },
    )
    
    // flye output channels
    renamed_contigs = renamed_contigs.mix(nanopore_contigs)
    assembly_log = assembly_log.mix(FLYE.out.log)
    raw_contigs = raw_contigs.mix(GUNZIP.out.gunzip)
    gfa = gfa.mix(FLYE.out.gfa)
    flye_info = flye_info.mix(FLYE.out.txt)
    ch_versions = ch_versions.mix(FLYE.out.versions)

    // illumina output channels
    renamed_contigs = renamed_contigs.mix(illumina_contigs)
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
    ch_flye_info       = flye_info
    versions           = ch_versions
}
