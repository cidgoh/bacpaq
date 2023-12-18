// INCLUDING MODULES
// Importing the necessary modules for the workflow
include { BAKTA_BAKTA } from '../../modules/nf-core/bakta/bakta/main'
include { BAKTA_BAKTADBDOWNLOAD } from '../../modules/nf-core/bakta/baktadbdownload/main'
include { PROKKA } from '../../modules/nf-core/prokka/main'

// Defining the workflow
workflow GENE_ANNOTATION{
    // Defining the input channel
    take:
        ch_genome

    // Defining the main process
    main:
        // Initializing empty channels for the output files
        prokka_gff = Channel.empty()
        bakta_gff = Channel.empty()
        ch_versions = Channel.empty()

        // Checking if the prodigal training file is provided
        if (!params.prodigal_training_file){
            prodigal = []
        }
        else{
            prodigals = file(params.prodigal_training_file, checkIfExists: true)
        }

        // Checking if the annotation protein file is provided
        if (!params.annotation_protein_file){
            proteins = []
        }
        else{
            proteins = file(params.annotation_protein_file, checkIfExists: true)
        }

        // Running PROKKA if skip_prokka is not set to true
        if (!params.skip_prokka){
            PROKKA(ch_genome, proteins, prodigal)
            prokka_gff = PROKKA.out.gff
            prokka_fna = PROKKA.out.fna
            prokka_faa = PROKKA.out.faa
            ch_versions   = ch_versions.mix(PROKKA.out.versions)
        }

        // Running BAKTA if skip_bakta is not set to true

        if(!params.skip_bakta){
            // Checking if BAKTA download is skipped
            if (!params.bakta_db){
                BAKTA_BAKTADBDOWNLOAD()
                bakta_db = BAKTA_BAKTADBDOWNLOAD.out.db
                ch_versions = ch_versions.mix(BAKTA_BAKTADBDOWNLOAD.out.versions)
            }
            else{
                bakta_db = Channel.from(params.bakta_db)

            }

            BAKTA_BAKTA(ch_genome, bakta_db, proteins, prodigal)
            bakta_gff = BAKTA_BAKTA.out.gff
            bakta_fna = BAKTA_BAKTA.out.fna
            bakta_faa = BAKTA_BAKTA.out.faa
            ch_versions = ch_versions.mix(BAKTA_BAKTA.out.versions)
        }

    // Defining the output channels
    emit:
        bakta_gff
        bakta_fna
        bakta_faa
        prokka_gff
        prokka_fna
        prokka_faa
        ch_versions
}
