include { ROARY } from '../../modules/nf-core/roary/main'

workflow PANGENOME_ANALYSIS {
    ch_gff = [
        [ id:'GCA_010673125.1', single_end:false ],
        [
            "/home/jyh25/scratch/hackathon/GCA_010673125.1.gff",
            "/home/jyh25/scratch/hackathon/GCA_010798115.1.gff"
        ]
    ]
    
    take:
    // should be ch_gff from prokka output
    // ch_gff

    main:
    ROARY(ch_gff)

    emit:
    results = ROARY.out.results
    versions = ROARY.out.versions
}