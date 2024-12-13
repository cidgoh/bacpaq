// import modules
include { CHECKM_LINEAGEWF  } from '../../modules/nf-core/checkm/lineagewf'
include { BUSCO             } from '../../modules/nf-core/busco'
include { QUAST             } from '../../modules/nf-core/quast'
include { CHECKLENGTH       } from '../../modules/local/checklength'

workflow ASSEMBLY_QC {
    take:
    assembly                        // ( val(meta), path(contigs.fasta) )

    main:
    ch_versions = Channel.empty()

    if (params.check_assembly_length) {
        CHECKLENGTH(assembly)
        ch_assembly_filtered = CHECKLENGTH.out.fasta_filtered
    }else{
        ch_assembly_filtered = assembly
    }

    if (!params.skip_checkm) {
    // RUN CHECKM
        ch_assembly_filtered
        .map { file(it[1]).getExtension() }
        .set { assembly_ext }

        CHECKM_LINEAGEWF(
            ch_assembly_filtered,
            assembly_ext,
            params.checkm_db != "null" ? params.checkm_db : assembly.map{ [] }
        )
        checkm_output = CHECKM_LINEAGEWF.out.checkm_output
        marker_file = CHECKM_LINEAGEWF.out.marker_file
        checkm_tsv = CHECKM_LINEAGEWF.out.checkm_tsv
        ch_versions = ch_versions.mix(CHECKM_LINEAGEWF.out.versions)
    }else{
        //empty output
        checkm_output = []
        marker_file = []
        checkm_tsv = []
    }

    // RUN QUAST
    if (!params.skip_quast) {
        if ( params.combine_quast ) {
            QUAST(
                ch_assembly_filtered.map{it[1]}.collect(),
                params.reference_genome_fasta != "null" ? params.reference_genome_fasta : [],
                params.reference_genome_gff != "null" ? params.reference_genome_gff : [],
                params.reference_genome_fasta != "null" ? ch_assembly_filtered.map{ true } : ch_assembly_filtered.map{ false },
                params.reference_genome_gff != "null" ? ch_assembly_filtered.map{ true } : ch_assembly_filtered.map{ false }
            )
        } else {
            QUAST(
                ch_assembly_filtered.map{it[1]},
                params.reference_genome_fasta != "null" ? params.reference_genome_fasta : ch_assembly_filtered.map{ [] },
                params.reference_genome_gff != "null" ? params.reference_genome_gff : ch_assembly_filtered.map{ [] },
                params.reference_genome_fasta != "null" ? ch_assembly_filtered.map{ true } : ch_assembly_filtered.map{ false },
                params.reference_genome_gff != "null" ? ch_assembly_filtered.map{ true } : ch_assembly_filtered.map{ false }
            )
        }
        ch_versions = ch_versions.mix(QUAST.out.versions)
        quast_results = QUAST.out.results
        quast_tsv = QUAST.out.tsv
    }else{
        //empty output
        quast_results = []
        quast_tsv = []
    }

    // RUN BUSCO
    if (!params.skip_busco) {
        BUSCO(
            ch_assembly_filtered,
            params.busco_lineage,
            params.busco_lineages_path != "null" ? params.busco_lineages_path : ch_assembly_filtered.map{ [] },
            params.busco_config ? params.busco_config : ch_assembly_filtered.map{ [] }
        )
        busco_batch_summary = BUSCO.out.batch_summary
        busco_short_summaries_txt = BUSCO.out.short_summaries_txt
        busco_short_summaries_json = BUSCO.out.short_summaries_json
        busco_dir = BUSCO.out.busco_dir
        ch_versions = ch_versions.mix(BUSCO.out.versions)
    }else{
        busco_batch_summary = []
        busco_short_summaries_txt = []
        busco_short_summaries_json = []
        busco_dir = []
    }

    emit:
    // CHECKM OUTPUTS
    checkm_output = checkm_output
    marker_file = marker_file
    checkm_tsv = checkm_tsv

    // QUAST OUTPUTS

    quast_results = quast_results
    quast_tsv = quast_tsv
    // BUSCO OUTPUTS
    busco_batch_summary = busco_batch_summary
    busco_short_summaries_txt = busco_short_summaries_txt
    busco_short_summaries_json = busco_short_summaries_json
    busco_dir = busco_dir
    versions = ch_versions


}
