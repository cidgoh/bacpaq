include { SNIPPY_RUN        } from '../../modules/nf-core/snippy/run'
include { SNIPPY_CORE       } from '../../modules/nf-core/snippy/core'
include { NUCMER            } from '../../modules/nf-core/nucmer'
include { GUBBINS           } from '../../modules/nf-core/gubbins'
include { MEDAKA_VARIANT    } from '../../modules/local/medaka/variant'
include { DELTA2VCF         } from '../../modules/local/delta2vcf'
include { CORESNPFILTER     } from '../../modules/local/coresnpfilter'
include { NUCMER_SAM        } from '../../modules/local/nucmersam'
include { SAMTOOLS_VIEW     } from '../../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX    } from '../../modules/nf-core/samtools/index'
include { SAMTOOLS_SORT     } from '../../modules/nf-core/samtools/sort'


workflow VARIANT_CALLING {

    take:
        ch_input

    main:
    // initialize channels
    ch_versions        = Channel.empty()
    ch_vcf_nucmer      = Channel.empty()
    ch_vcf_snippy      = Channel.empty()
    ch_vcf_medaka      = Channel.empty()
    ch_vcf_gubbins     = Channel.empty()
    ch_sam_nucmer      = Channel.empty()
    ch_bam_nucmer      = Channel.empty()
    ch_bam_snippy      = Channel.empty()
    ch_bam_medaka      = Channel.empty()
    ch_bai_snippy      = Channel.empty()
    ch_bai_medaka      = Channel.empty()
    ch_bai_nucmer      = Channel.empty()
    ch_full_aln        = Channel.empty()
    ch_core_aln        = Channel.empty()

    // branch input channel according to meta.mode
    ch_input
        .branch { meta, path ->
            illumina: meta.mode == 'illumina'
            nanopore: meta.mode == 'nanopore'
            genome: true
        }
        .set { ch_input_seq }
    // reference genome
    Channel
        .fromPath(params.reference_genome)
        .set { ch_reference_genome }

    // GENOME: RUN NUCMER
    ch_input_seq.genome
        .combine(ch_reference_genome)
        .map { meta, genome, ref -> [ meta, ref, genome ] }
        .set { ch_nucmer }

    NUCMER(ch_nucmer)
    ch_versions    = ch_versions.mix(NUCMER.out.versions)
    NUCMER_SAM(ch_nucmer)
    ch_versions    = ch_versions.mix(NUCMER_SAM.out.versions)  
    ch_sam_nucmer  = NUCMER_SAM.out.sam_fixed
    SAMTOOLS_VIEW(
        ch_sam_nucmer.map { meta, sam -> [ meta, sam, [] ] },
        ch_reference_genome.first(), 
        []
    )
    ch_bam_nucmer  = SAMTOOLS_VIEW.out.bam
    ch_versions    = ch_versions.mix(SAMTOOLS_VIEW.out.versions)
    SAMTOOLS_SORT(ch_bam_nucmer)
    ch_sbam_nucmer = SAMTOOLS_SORT.out.bam
    ch_versions    = ch_versions.mix(SAMTOOLS_SORT.out.versions)
    SAMTOOLS_INDEX(ch_sbam_nucmer)
    ch_bai_nucmer  = SAMTOOLS_INDEX.out.bai
    ch_versions    = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
    DELTA2VCF(NUCMER.out.delta)
    ch_versions    = ch_versions.mix(DELTA2VCF.out.versions)
    ch_vcf_nucmer  = DELTA2VCF.out.vcf
    

    // NANOPORE: RUN MEDAKA
    ch_medaka  = ch_input_seq.nanopore
    MEDAKA_VARIANT(ch_medaka, ch_reference_genome.first())
    ch_vcf_medaka = MEDAKA_VARIANT.out.alignment_vcf
    ch_bam_medaka = MEDAKA_VARIANT.out.bam
    ch_bai_medaka = MEDAKA_VARIANT.out.bam_bai

    // ILLUMINA: RUN SNIPPY
    ch_snippy = ch_input_seq.illumina
    SNIPPY_RUN(ch_snippy, ch_reference_genome.first())
    ch_vcf_snippy   = SNIPPY_RUN.out.vcf
    ch_bam_snippy   = SNIPPY_RUN.out.bam
    ch_bai_snippy   = SNIPPY_RUN.out.bai
    
    ch_snippy_aligned_fa = SNIPPY_RUN.out.aligned_fa
        .map { it[1] }
        .collect()
        .map { [[ id: 'core_aln' ], it ] }

    ch_snippy_vcf=SNIPPY_RUN.out.vcf
        .map { it[1] }
        .collect()
        .map { [ [id: 'core_aln'], it ] }

    // generate core SNP alignment
    ch_snippy_core = ch_snippy_vcf.combine(ch_snippy_aligned_fa, by: 0)    
    SNIPPY_CORE(ch_snippy_core, ch_reference_genome)
    ch_versions = ch_versions.mix(SNIPPY_CORE.out.versions)
    ch_full_aln = SNIPPY_CORE.out.full_aln

    // filter recombinant sites
    if ( !params.skip_gubbins) {

        ch_gubbins = SNIPPY_CORE.out.full_aln
            .map { it[1] }
        GUBBINS(ch_gubbins)
        ch_vcf_gubbins      = GUBBINS.out.vcf
        ch_versions = ch_versions.mix(GUBBINS.out.versions)
        ch_coresnpfilter = GUBBINS.out.fasta

    } else {

        ch_coresnpfilter = SNIPPY_CORE.out.full_aln
        ch_core_aln      = SNIPPY_CORE.out.aln

    }

    // filter SNP alignment by core SNP thresholds (e.g., 0.99)
    if ( !params.skip_coresnpfilter ) {

        CORESNPFILTER(ch_coresnpfilter)
        ch_versions = ch_versions.mix(CORESNPFILTER.out.versions)
        ch_core_aln = CORESNPFILTER.out.alignment

    }

    emit:
    versions            = ch_versions
    vcf_nucmer          = ch_vcf_nucmer
    vcf_snippy          = ch_vcf_snippy
    vcf_medaka          = ch_vcf_medaka
    vcf_gubbins         = ch_vcf_gubbins
    sam_nucmer          = ch_sam_nucmer
    bam_snippy          = ch_bam_snippy
    bam_medaka          = ch_bam_medaka
    bam_nucmer          = ch_bam_nucmer 
    bam_nucmer_sorted   = ch_sbam_nucmer
    bai_snippy          = ch_bai_snippy
    bai_medaka          = ch_bai_medaka
    bai_nucmer          = ch_bai_nucmer
    full_aln            = ch_full_aln
    core_aln            = ch_core_aln
}
