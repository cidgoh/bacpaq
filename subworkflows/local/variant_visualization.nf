include { IGVREPORTS } from '../../modules/nf-core/igvreports/main'
include { IQTREE } from '../../modules/nf-core/iqtree/main'

workflow VARIANT_VISUALIZATION{


    ch_vcf = Channel
        .fromPath("/home/myee/projects/def-sponsor00/myee/bacpaq/null/snippy/SRR1090340{1,2}/SRR1090340{1,2}.vcf")
        .collect()
        .map { [[ id: 'igv' ], it ] }

    ch_bam=Channel
        .fromPath("/home/myee/projects/def-sponsor00/myee/bacpaq/null/snippy/SRR1090340*/SRR1090340*.bam")
        .collect()
        .map { [[ id: 'igv' ], it ] }

    ch_bam_index=Channel
        .fromPath("/home/myee/projects/def-sponsor00/myee/bacpaq/null/snippy/SRR1090340*/SRR1090340*.bam.bai")
        .collect()
        .map { [[ id: 'igv' ], it ] }
    
    ch_igv_input_1=ch_vcf.combine(ch_bam, by: 0)
    
    ch_igv_input_2=ch_igv_input_1.combine(ch_bam_index, by: 0)

    ch_reference = Channel.fromPath("/mnt/cidgoh-object-storage/hackathon/variant_calls/reference/GCF_009858895.2.fasta")
        .map { [[ id: 'reference' ], it ] }

    ch_reference_index=Channel.fromPath("/mnt/cidgoh-object-storage/hackathon/variant_calls/reference/GCF_009858895.2.fasta.fai")
            .map { [[ id: 'reference' ], it ] }
    
    ch_ref_and_index=ch_reference.combine(ch_reference_index, by:0)

    //IGVREPORTS(ch_igv_input_2, ch_ref_and_index)

    ch_alignment=Channel.fromPath("/mnt/cidgoh-object-storage/hackathon/variant_calls/clean.core.aln")
        .map { [[id:"tree"], it,[]] }

    IQTREE(ch_alignment,[],[],[],[],[],[],[],[],[],[],[],[])


    



}