//
// VARIANT DETECTION WORKFLOW:
// Accepts short, long, assembled sequence data to identify mutations
// differing from reference alleles based on a user-supplied reference
// genome. Visualization of identified polymorphisms is supported.
//

// import subworkflows
include { VARIANT_CALLING } from './variantcalling'
include { VARIANT_VIS     } from './variantvis'

workflow VARIANT_DETECTION {
    take:
    reads  // channel containing path to sequence data in FASTQ format
    genome // channel containing path to reference genome in FASTA format

    main:
    // initialize channels
    ch_versions = Channel.empty()
    ch_vcf = Channel.empty()
    ch_bam = Channel.empty()
    ch_bai = Channel.empty()
    ch_aln_fa = Channel.empty()
    ch_igvreports_vcf = Channel.empty()
    ch_igvreports_bam = Channel.empty()
    ch_iqtree_nwk = Channel.empty()
    ch_iqtree_report = Channel.empty()

    // VARIANT CALLING SUBWORKFLOW
    VARIANT_CALLING(reads, genome)
    ch_versions = ch_versions.mix(VARIANT_CALLING.out.versions)

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
            ch_aln_fa,
        )

        ch_versions = ch_versions.mix(VARIANT_VIS.out.versions)
    }

    emit:
    versions      = ch_versions // channel: [ versions.yml ]
    vcf           = ch_vcf // channel: [ val(meta), [ *.vcf ] ] from snippy, nucmer, or medaka
    core_aln      = ch_aln_fa // channel: [ val(meta), [ *.aln ] ] from core-snp-filter
    bam           = ch_bam // channel: [ val(meta), [ *.bam ] ]
    bai           = ch_bai // channel: [ val(meta), [ *.bai ] ]
    igv_vcf       = ch_igvreports_vcf // channel: [ val(meta), [ *.html ] ]
    igv_bam       = ch_igvreports_bam // channel: [ val(meta), [ *.html ] ]
    iqtree_nwk    = ch_iqtree_nwk // channel: [ val(meta), [ *.treefile] ]
    iqtree_report = ch_iqtree_report // channel: [ val(meta), [ *.iqtree] ]
}
