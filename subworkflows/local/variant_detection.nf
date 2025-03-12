// 
// VARIANT DETECTION WORKFLOW:
// Accepts short, long, assembled sequence data to identify mutations
// differing from reference alleles based on a user-supplied reference
// genome. Visualization of identified polymorphisms is supported.
//

// import subworkflows
include { VARIANT_CALLING } from './variantcalling'
include { VARIANT_VIS } from './variantvis'

workflow VARIANT_DETECTION {

    take: ch_input // channel containing path to sequence data in FASTQ or FASTA format
    
    main:
        
        // VARIANT CALLING SUBWORKFLOW
        VARIANT_CALLING(ch_input)        
        
        // VARIANT VIZ SUBWORKFLOW
        if (!params.skip_variant_viz) {

            ch_vcf = VARIANT_CALLING.out.vcf_snippy
                .concat(VARIANT_CALLING.out.vcf_medaka)
                .concat(VARIANT_CALLING.out.vcf_nucmer)
            ch_bam = VARIANT_CALLING.out.bam_snippy
                .concat(VARIANT_CALLING.out.bam_medaka)
                .concat(VARIANT_CALLING.out.bam_nucmer_sorted)
            ch_bai = VARIANT_CALLING.out.bai_snippy
                .concat(VARIANT_CALLING.out.bai_medaka)
                .concat(VARIANT_CALLING.out.bai_nucmer)
            ch_aln_fa = VARIANT_CALLING.out.core_aln

            VARIANT_VIS(
                ch_vcf,
                ch_bam,
                ch_bai,
                ch_aln_fa
            )
        }

}
