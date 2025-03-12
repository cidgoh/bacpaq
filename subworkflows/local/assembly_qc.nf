// import modules
include { CHECKM_LINEAGEWF } from '../../modules/nf-core/checkm/lineagewf'
include { BUSCO            } from '../../modules/nf-core/busco'
include { QUAST            } from '../../modules/nf-core/quast'

workflow ASSEMBLY_QC {
    take:
    assembly

    main:
    ch_versions = Channel.empty()


    if (!params.skip_checkm) {
        // RUN CHECKM
        assembly
            .map { file(it[1]).getExtension() }
            .set { assembly_ext }

        CHECKM_LINEAGEWF(
            assembly,
            assembly_ext,
            params.checkm_db != "null" ? params.checkm_db : assembly.map { [] }
        )
        checkm_output = CHECKM_LINEAGEWF.out.checkm_output
        marker_file = CHECKM_LINEAGEWF.out.marker_file
        checkm_tsv = CHECKM_LINEAGEWF.out.checkm_tsv
        ch_versions = ch_versions.mix(CHECKM_LINEAGEWF.out.versions)
    }
    else {
        //empty output
        checkm_output = []
        marker_file = []
        checkm_tsv = []
    }

    // RUN QUAST
    if (!params.skip_quast) {
        if (params.combine_quast) {
            combined_assemblies = assembly
                .flatten()
                .collect()
                .map { assemblies ->
                    def meta = [id: 'combined_assemblies']
                    def all_files = assemblies.collect { it[1] }.flatten()
                    tuple(meta, all_files)
                }
            QUAST(
                combined_assemblies,
                params.reference_genome_fasta != "null" ? params.reference_genome_fasta : [[], []],
                params.reference_genome_gff != "null" ? params.reference_genome_gff : [[], []]
            )
        }
        else {
            QUAST(
                assembly,
                params.reference_genome_fasta != "null" ? params.reference_genome_fasta : [[], []],
                params.reference_genome_gff != "null" ? params.reference_genome_gff : [[], []]
            )
        }
        ch_versions = ch_versions.mix(QUAST.out.versions)
        quast_results = QUAST.out.results
        quast_tsv = QUAST.out.tsv
    }
    else {
        //empty output
        quast_results = []
        quast_tsv = []
    }

    // RUN BUSCO
    if (!params.skip_busco) {
        BUSCO(
            assembly,
            params.busco_lineage,
            params.busco_lineages_path != "null" ? params.busco_lineages_path : assembly.map { [] },
            params.busco_config ? params.busco_config : assembly.map { [] }
        )
        busco_batch_summary = BUSCO.out.batch_summary
        busco_short_summaries_txt = BUSCO.out.short_summaries_txt
        busco_short_summaries_json = BUSCO.out.short_summaries_json
        busco_dir = BUSCO.out.busco_dir
        ch_versions = ch_versions.mix(BUSCO.out.versions)
    }
    else {
        busco_batch_summary = []
        busco_short_summaries_txt = []
        busco_short_summaries_json = []
        busco_dir = []
    }

    emit:
    checkm_output              = checkm_output
    marker_file                = marker_file
    checkm_tsv                 = checkm_tsv
    quast_results              = quast_results
    quast_tsv                  = quast_tsv
    busco_batch_summary        = busco_batch_summary
    busco_short_summaries_txt  = busco_short_summaries_txt
    busco_short_summaries_json = busco_short_summaries_json
    busco_dir                  = busco_dir
    versions                   = ch_versions
}
