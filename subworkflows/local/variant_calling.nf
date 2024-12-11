// import modules

include { SNIPPY_CORE } from '../../modules/nf-core/snippy/core/main'
include { SNIPPY_RUN } from '../../modules/nf-core/snippy/run/main'
include { GUBBINS } from '../../modules/nf-core/gubbins/main'


workflow VARIANT_CALLING{
        
   
    ch_reference = Channel.fromPath("/mnt/cidgoh-object-storage/hackathon/variant_calls/reference/GCF_009858895.2.fasta")
    ch_reads = Channel
    .fromPath("/mnt/cidgoh-object-storage/hackathon/variant_calls/short_reads/SRR1090340*/SRR1090340*_{1,2}.fastq.gz")
    .map { file -> 
        def meta = [id: file.simpleName.replaceAll(/_[12]$/, '')]
        return tuple(meta, file)
    }
    .groupTuple(by: 0)


    
    SNIPPY_RUN(ch_reads,ch_reference.first())

    ch_aligned_fa=SNIPPY_RUN.out.aligned_fa
        .map { it[1] }
        .collect()
        .map { [[ id: 'core_aln' ], it ] }

    ch_vcf=SNIPPY_RUN.out.vcf
        .map { it[1] }
        .collect()
        .map { [ [id: 'core_aln'], it ] }


    
    ch_core_input=ch_vcf.combine(ch_aligned_fa, by: 0)
    
    
    SNIPPY_CORE(ch_core_input, ch_reference)

    ch_gubbins=SNIPPY_CORE.out.full_aln
        .map { it[1 ]}

    GUBBINS(ch_gubbins)

}