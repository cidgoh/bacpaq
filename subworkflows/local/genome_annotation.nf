/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
//WorkflowSeqqc.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
//def checkPathParamList = [  params.multiqc_config, params.fasta ]
//for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
// if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
//ch_input = file(params.input)
/*
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
//include { INPUT_CHECK           } from './input_check'
include { GENE_ANNOTATION    } from './gene_annotation'
include { PHAGE              } from './phage'
include { PANGENOME_ANALYSIS } from './pangenome_analysis'
include { PLASMIDS           } from './plasmids'
include { AMR_ANNOTATION     } from './amr_annotation'
include { CRISPRS            } from './crisprs'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
//include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



workflow ANNOTATION {
    take:
    ch_genome

    main:
    ch_versions = Channel.empty()
    ch_pangenome_gff = Channel.empty()
    ch_multiqc_files = Channel.empty()

    if (!params.skip_amr_annotation) {
        // Annotate AMR using ABRICATE
        AMR_ANNOTATION(ch_genome)
        ch_versions = ch_versions.mix(AMR_ANNOTATION.out.versions)
    }

    //phage_input = INPUT_CHECK.out.reads
    if (!params.skip_phage_annotation) {
        // Annotate phages using VIRSORTER2
        PHAGE(ch_genome)
        ch_versions = ch_versions.mix(PHAGE.out.versions)
    }

    //plasmid_input = INPUT_CHECK.out.reads
    if (!params.skip_plasmid_analysis) {
        // Run plasmid annotation
        PLASMIDS(ch_genome)
        ch_versions = ch_versions.mix(PLASMIDS.out.versions)
    }

    //crispr_input = INPUT_CHECK.out.reads
    if (!params.skip_crispr_analysis) {
        // Run CRISPR analysis
        CRISPRS(ch_genome)
        ch_versions = ch_versions.mix(CRISPRS.out.versions)
    }

    //gene_annotation_input = INPUT_CHECK.out.reads
    if (!params.skip_gene_annotation) {
        // Annotate genomes using Prokka and bakta
        GENE_ANNOTATION(ch_genome)
        ch_multiqc_files = ch_multiqc_files.mix(GENE_ANNOTATION.out.prokka_txt.collect { it[1] }.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(GENE_ANNOTATION.out.bakta_txt.collect { it[1] }.ifEmpty([]))
        ch_versions = ch_versions.mix(GENE_ANNOTATION.out.versions)

        if (!params.skip_pangenome_analysis) {
            // Run pangenome analysis using bakta output

            if (params.pangenome_input == "bakta") {
                gffs = GENE_ANNOTATION.out.bakta_gff
                gffs
                    .map { [it[1]] }
                    .collect()
                    .map { gff -> [[id: "bakta_pangenome"], gff] }
                    .set { ch_pangenome_gff }
            }
            else {
                gffs = GENE_ANNOTATION.out.prokka_gff
                gffs
                    .map { [it[1]] }
                    .collect()
                    .map { gff -> [[id: "prokka_pangenome"], gff] }
                    .set { ch_pangenome_gff }
            }

            PANGENOME_ANALYSIS(ch_pangenome_gff)
            ch_versions = ch_versions.mix(PANGENOME_ANALYSIS.out.versions)
        }
    }

    emit:
    versions      = ch_versions
    multiqc_files = ch_multiqc_files
}
