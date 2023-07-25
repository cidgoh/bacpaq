// import modules
include { CHECKM_LINEAGEWF } from '../../modules/nf-core/checkm/lineagewf'
include { BUSCO } from '../../modules/nf-core/busco'
include { QUAST } from '../../modules/nf-core/quast'

workflow ASSEMBLY_QC {
    take: 
    assembly

    main:
    if (!params.skip_checkm) {
    // RUN CHECKM
        assembly
        .map { file(it[1]).getExtension() }
        .set { assembly_ext }

        CHECKM_LINEAGEWF(
            assembly,
            assembly_ext,   
            params.checkm_db ? params.checkm_db : []
        )
        checkm_output = CHECKM_LINEAGEWF.out.checkm_output
        marker_file = CHECKM_LINEAGEWF.out.marker_file
        checkm_tsv = CHECKM_LINEAGEWF.out.checkm_tsv
        checkm_versions = CHECKM_LINEAGEWF.out.versions    
    }else{
        //empty output
        checkm_output = []
        marker_file = []
        checkm_tsv = []
        checkm_versions = []    
    }
    
    // RUN QUAST
    if (!params.skip_quast) {
        if ( params.combine_quast ) {
            QUAST(
                assembly.map{it[1]}.collect(),
                params.reference_genome_fasta != "null" ? params.reference_genome_fasta : [],
                params.reference_genome_gff != "null" ? params.reference_genome_gff : [],
                params.reference_genome_fasta != "null" ? assembly.map{ true } : assembly.map{ false },
                params.reference_genome_gff != "null" ? assembly.map{ true } : assembly.map{ false }
            )
        } else {
            QUAST(
                assembly.map{it[1]},
                params.reference_genome_fasta != "null" ? params.reference_genome_fasta : assembly.map{ [] },
                params.reference_genome_gff != "null" ? params.reference_genome_gff : assembly.map{ [] },
                params.reference_genome_fasta != "null" ? assembly.map{ true } : assembly.map{ false },
                params.reference_genome_gff != "null" ? assembly.map{ true } : assembly.map{ false }
            )
        }
        quast_verions = QUAST.out.versions
        quast_results = QUAST.out.results
        quast_tsv = QUAST.out.tsv
    }else{
        //empty output
        quast_verions = []
        quast_results = []
        quast_tsv = []    
    }

    // RUN BUSCO
    if (!params.skip_busco) {
        BUSCO(
            assembly,
            params.busco_lineage,
            params.busco_lineages_path != "null" ? params.busco_lineages_path : assembly.map{ [] },
            params.busco_config ? params.busco_config : assembly.map{ [] }
        )
        busco_batch_summary = BUSCO.out.batch_summary
        busco_short_summaries_txt = BUSCO.out.short_summaries_txt 
        busco_short_summaries_json = BUSCO.out.short_summaries_json
        busco_dir = BUSCO.out.busco_dir
        busco_versions = BUSCO.out.versions
    }else{
        busco_batch_summary = []
        busco_short_summaries_txt = []
        busco_short_summaries_json = []
        busco_dir = []
        busco_versions = []
    }

    emit:
    // CHECKM OUTPUTS
    checkm_output = checkm_output
    marker_file = marker_file
    checkm_tsv = checkm_tsv
    checkm_versions = checkm_versions   
    // QUAST OUTPUTS
    quast_verions = quast_verions
    quast_results = quast_results
    quast_tsv = quast_tsv
    // BUSCO OUTPUTS
    busco_batch_summary = busco_batch_summary
    busco_short_summaries_txt = busco_short_summaries_txt 
    busco_short_summaries_json = busco_short_summaries_json
    busco_dir = busco_dir
    busco_versions = busco_versions

}