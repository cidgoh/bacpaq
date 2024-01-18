// INCLUDE MODULES
include { VIRSORTER2_RUN } from '../../modules/local/virsorter2/run/main.nf'  
include { VIRSORTER2_SETUP } from '../../modules/local/virsorter2/setup/main.nf'  

// Phage ID workflow
workflow PHAGE {
    take: genome
    main:
        
        // inititalize workflow output channels
        ch_versions = Channel.empty()
        ch_virsorter2_fa = Channel.empty()
        ch_virsorter2_boundary = Channel.empty()

        if (!params.skip_virsorter2) {
            // download virsorter2 db if db path is not provided
            // and running in conda envs
            if ( !params.virsorter_db & !workflow.containerEngine ) {
                date = new Date()
                date_ymd = date.format("yyyyMMdd")
                db_prefix = "virsorter2_db"
                db_name = [db_prefix, date_ymd].join('_')
                ch_db_meta = [id:db_name]
                VIRSORTER2_SETUP(ch_db_meta)
                
                ch_db_path = VIRSORTER2_SETUP.out.db
                ch_versions = ch_versions.mix(VIRSORTER2_SETUP.out.versions)
            } else {
                ch_db_path = []
            }
            
            // run phage prediction with virsorter2
            VIRSORTER2_RUN(genome, ch_db_path)
            ch_versions = ch_versions.mix(VIRSORTER2_RUN.out.versions)
            ch_virsorter2_fa = VIRSORTER2_RUN.out.fasta
            ch_virsorter2_boundary = VIRSORTER2_RUN.out.boundary
        }
                
    emit:
        versions = ch_versions
        virsorter2_fa = ch_virsorter2_fa
        virsorter2_boundary = ch_virsorter2_boundary
}