include { SNIPPY_RUN     } from '../../modules/nf-core/snippy/run'
include { SNIPPY_CORE    } from '../../modules/nf-core/snippy/core'
include { NUCMER         } from '../../modules/nf-core/nucmer'
include { GUBBINS        } from '../../modules/nf-core/gubbins'
include { MEDAKA_VARIANT } from '../../modules/local/medaka/variant'
include { DELTA2VCF      } from '../../modules/local/delta2vcf'
include { CORESNPFILTER  } from '../../modules/local/coresnpfilter'
include { NUCMER_SAM     } from '../../modules/local/nucmersam'
include { SAMTOOLS_VIEW  } from '../../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index'
include { SAMTOOLS_SORT  } from '../../modules/nf-core/samtools/sort'
include { TABIX_TABIX as TABIX_SNIPPY; TABIX_TABIX as TABIX_MEDAKA; TABIX_TABIX as TABIX_NUCMER; TABIX_TABIX as TABIX_CORE } from '../../modules/nf-core/tabix/tabix'
include { TABIX_BGZIP as BGZIP_SNIPPY; TABIX_BGZIP as BGZIP_MEDAKA; TABIX_BGZIP as BGZIP_NUCMER; TABIX_BGZIP as BGZIP_CORE } from '../../modules/nf-core/tabix/bgzip'


workflow VARIANT_CALLING {
    take:
    ch_reads  // channel containing path to sequence data in FASTQ format
    ch_genome // channel containing path to reference genome in FASTA format

    main:
    // initialize channels
    ch_versions = Channel.empty()
    ch_vcf_nucmer = Channel.empty()
    ch_vcf_snippy = Channel.empty()
    ch_vcf_medaka = Channel.empty()
    ch_vcf_gubbins = Channel.empty()
    ch_vcf_bgz_nucmer = Channel.empty()
    ch_vcf_bgz_snippy = Channel.empty()
    ch_vcf_bgz_medaka = Channel.empty()
    ch_vcf_gubbins = Channel.empty()
    ch_sam_nucmer = Channel.empty()
    ch_bam_nucmer = Channel.empty()
    ch_bam_snippy = Channel.empty()
    ch_bam_medaka = Channel.empty()
    ch_bai_snippy = Channel.empty()
    ch_bai_medaka = Channel.empty()
    ch_bai_nucmer = Channel.empty()
    ch_vci_snippy = Channel.empty()
    ch_vci_medaka = Channel.empty()
    ch_vci_nucmer = Channel.empty()
    // snippy output channels
    ch_tab_snippy = Channel.empty()
    ch_csv_snippy = Channel.empty()
    ch_html_snippy = Channel.empty()
    ch_bed_snippy = Channel.empty()
    ch_gff_snippy = Channel.empty()
    ch_log_snippy = Channel.empty()
    ch_aligned_fa_snippy = Channel.empty()
    ch_consensus_fa_snippy = Channel.empty()
    ch_consensus_subs_fa_snippy = Channel.empty()
    ch_raw_vcf_snippy = Channel.empty()
    ch_filt_vcf_snippy = Channel.empty()
    ch_vcf_csi_snippy = Channel.empty()
    ch_vcf_gz_snippy = Channel.empty()
    ch_vcf_csi_snippy = Channel.empty()
    ch_vcf_gz_snippy = Channel.empty()
    ch_txt_snippy = Channel.empty()
    // snippy core output channels
    ch_full_aln = Channel.empty()
    ch_core_aln = Channel.empty()
    ch_core_tab = Channel.empty()
    ch_core_vcf = Channel.empty()
    ch_core_vcf_bgz = Channel.empty()
    ch_core_txt = Channel.empty()
    ch_core_vci = Channel.empty()

    // branch input channel according to meta.mode
    ch_reads
        .branch { meta, path ->
            illumina: meta.mode == 'illumina'
            nanopore: meta.mode == 'nanopore'
        }
        .set { ch_input_seq }
    // reference genome
    Channel.fromPath(params.reference_genome)
        .set { ch_reference_genome }

    // GENOME: RUN NUCMER
    ch_genome
        .combine(ch_reference_genome)
        .map { meta, genome, ref -> [meta, ref, genome] }
        .set { ch_nucmer }
    NUCMER(ch_nucmer)
    ch_versions = ch_versions.mix(NUCMER.out.versions)
    NUCMER_SAM(ch_nucmer)
    ch_versions = ch_versions.mix(NUCMER_SAM.out.versions)
    ch_sam_nucmer = NUCMER_SAM.out.sam_fixed
    reference_genome = ch_reference_genome.map { path ->
        def meta = path.getName().replaceFirst(/\.[^.]+$/, '')
        // Extract filename without extension
        tuple(meta, path)
    }
    SAMTOOLS_VIEW(
        ch_sam_nucmer.map { meta, sam -> [meta, sam, []] },
        reference_genome,
        [],
    )
    ch_bam_nucmer = SAMTOOLS_VIEW.out.bam
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)
    SAMTOOLS_SORT(ch_bam_nucmer, reference_genome)
    ch_sbam_nucmer = SAMTOOLS_SORT.out.bam
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)
    SAMTOOLS_INDEX(ch_sbam_nucmer)
    ch_bai_nucmer = SAMTOOLS_INDEX.out.bai
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
    DELTA2VCF(NUCMER.out.delta)
    ch_versions = ch_versions.mix(DELTA2VCF.out.versions)
    ch_vcf_nucmer = DELTA2VCF.out.vcf
    BGZIP_NUCMER(ch_vcf_nucmer) // compress the VCF file
    ch_versions = ch_versions.mix(BGZIP_NUCMER.out.versions)
    ch_vcf_bgz_nucmer = BGZIP_NUCMER.out.output
    TABIX_NUCMER(ch_vcf_bgz_nucmer) // index the VCF file
    ch_versions = ch_versions.mix(TABIX_NUCMER.out.versions)
    ch_vci_nucmer = TABIX_NUCMER.out.tbi 

    // NANOPORE: RUN MEDAKA
    ch_medaka = ch_input_seq.nanopore
    MEDAKA_VARIANT(ch_medaka, ch_reference_genome.first())
    ch_vcf_medaka = MEDAKA_VARIANT.out.alignment_vcf
    ch_bam_medaka = MEDAKA_VARIANT.out.bam
    ch_bai_medaka = MEDAKA_VARIANT.out.bam_bai
    BGZIP_MEDAKA(ch_vcf_medaka) // compress the VCF file
    ch_versions = ch_versions.mix(BGZIP_MEDAKA.out.versions)
    ch_vcf_bgz_medaka = BGZIP_MEDAKA.out.output
    TABIX_MEDAKA(ch_vcf_bgz_medaka) // index the VCF file
    ch_versions = ch_versions.mix(TABIX_MEDAKA.out.versions)
    ch_vci_medaka = TABIX_MEDAKA.out.tbi

    // ILLUMINA: RUN SNIPPY
    ch_snippy = ch_input_seq.illumina
    SNIPPY_RUN(ch_snippy, ch_reference_genome.first())
    ch_vcf_snippy = SNIPPY_RUN.out.vcf
    ch_bam_snippy = SNIPPY_RUN.out.bam
    ch_bai_snippy = SNIPPY_RUN.out.bai
    ch_tab_snippy = SNIPPY_RUN.out.tab
    ch_csv_snippy = SNIPPY_RUN.out.csv
    ch_html_snippy = SNIPPY_RUN.out.html
    ch_bed_snippy = SNIPPY_RUN.out.bed
    ch_gff_snippy = SNIPPY_RUN.out.gff
    ch_log_snippy = SNIPPY_RUN.out.log
    ch_aligned_fa_snippy = SNIPPY_RUN.out.aligned_fa
    ch_consensus_fa_snippy = SNIPPY_RUN.out.consensus_fa
    ch_consensus_subs_fa_snippy = SNIPPY_RUN.out.consensus_subs_fa
    ch_raw_vcf_snippy = SNIPPY_RUN.out.raw_vcf
    ch_filt_vcf_snippy = SNIPPY_RUN.out.filt_vcf
    ch_vcf_csi_snippy = SNIPPY_RUN.out.vcf_csi
    ch_vcf_gz_snippy = SNIPPY_RUN.out.vcf_gz
    ch_vcf_csi_snippy = SNIPPY_RUN.out.vcf_csi
    ch_vcf_gz_snippy = SNIPPY_RUN.out.vcf_gz
    ch_txt_snippy = SNIPPY_RUN.out.txt

    BGZIP_SNIPPY(ch_vcf_snippy) // compress the VCF file
    ch_versions = ch_versions.mix(BGZIP_SNIPPY.out.versions)
    ch_vcf_bgz_snippy = BGZIP_SNIPPY.out.output
    TABIX_SNIPPY(ch_vcf_bgz_snippy) // index the VCF file
    ch_versions = ch_versions.mix(TABIX_SNIPPY.out.versions)
    ch_vci_snippy = TABIX_SNIPPY.out.tbi

    ch_snippy_aligned_fa = SNIPPY_RUN.out.aligned_fa
        .map { it[1] }
        .collect()
        .map { [[id: 'core_aln'], it] }

    ch_snippy_vcf = SNIPPY_RUN.out.vcf
        .map { it[1] }
        .collect()
        .map { [[id: 'core_aln'], it] }

    // generate core SNP alignment
    ch_snippy_core = ch_snippy_vcf.combine(ch_snippy_aligned_fa, by: 0)
    SNIPPY_CORE(ch_snippy_core, ch_reference_genome)
    ch_versions = ch_versions.mix(SNIPPY_CORE.out.versions)
    ch_full_aln = SNIPPY_CORE.out.full_aln
    ch_core_tab = SNIPPY_CORE.out.tab
    ch_core_vcf = SNIPPY_CORE.out.vcf
    ch_core_txt = SNIPPY_CORE.out.txt

    BGZIP_CORE(ch_core_vcf) // compress the VCF file
    ch_versions = ch_versions.mix(BGZIP_CORE.out.versions)
    ch_core_vcf_bgz = BGZIP_CORE.out.output
    TABIX_CORE(ch_core_vcf_bgz) // index the VCF file
    ch_versions = ch_versions.mix(TABIX_CORE.out.versions)
    ch_core_vci = TABIX_CORE.out.tbi

    // filter recombinant sites
    if (!params.skip_gubbins) {
        ch_gubbins = SNIPPY_CORE.out.full_aln.map { it[1] }
        GUBBINS(ch_gubbins)
        ch_vcf_gubbins = GUBBINS.out.vcf
        ch_versions = ch_versions.mix(GUBBINS.out.versions)
        ch_coresnpfilter = GUBBINS.out.fasta.map { fasta -> tuple([id: 'core_aln'], fasta) }
    }
    else {
        ch_coresnpfilter = SNIPPY_CORE.out.full_aln
        ch_core_aln = SNIPPY_CORE.out.aln
    }

    // filter SNP alignment by core SNP thresholds (e.g., 0.99)
    if (!params.skip_coresnpfilter) {
        CORESNPFILTER(ch_coresnpfilter)
        ch_versions = ch_versions.mix(CORESNPFILTER.out.versions)
        ch_core_aln = CORESNPFILTER.out.alignment
    }

    emit:
    versions          = ch_versions
    vcf_nucmer        = ch_vcf_nucmer
    vcf_snippy        = ch_vcf_snippy
    vcf_medaka        = ch_vcf_medaka
    vcf_gubbins       = ch_vcf_gubbins
    vcf_bgz_nucmer        = ch_vcf_bgz_nucmer
    vcf_bgz_snippy        = ch_vcf_bgz_snippy
    vcf_bgz_medaka        = ch_vcf_bgz_medaka
    sam_nucmer        = ch_sam_nucmer
    bam_snippy        = ch_bam_snippy
    bam_medaka        = ch_bam_medaka
    bam_nucmer        = ch_bam_nucmer
    bam_nucmer_sorted = ch_sbam_nucmer
    bai_snippy        = ch_bai_snippy
    bai_medaka        = ch_bai_medaka
    bai_nucmer        = ch_bai_nucmer
    vci_snippy        = ch_vci_snippy
    vci_medaka        = ch_vci_medaka
    vci_nucmer        = ch_vci_nucmer
    // snippy outputs
    tab_snippy        = ch_tab_snippy
    csv_snippy        = ch_csv_snippy
    html_snippy       = ch_html_snippy
    bed_snippy        = ch_bed_snippy
    gff_snippy        = ch_gff_snippy
    log_snippy        = ch_log_snippy
    aligned_fa_snippy = ch_aligned_fa_snippy
    consensus_fa_snippy = ch_consensus_fa_snippy
    consensus_subs_fa_snippy = ch_consensus_subs_fa_snippy
    raw_vcf_snippy    = ch_raw_vcf_snippy
    filt_vcf_snippy   = ch_filt_vcf_snippy
    vcf_csi_snippy    = ch_vcf_csi_snippy
    vcf_gz_snippy     = ch_vcf_gz_snippy
    txt_snippy        = ch_txt_snippy
    // snippy core outputs
    core_tab          = ch_core_tab
    core_vcf          = ch_core_vcf
    core_txt          = ch_core_txt
    full_aln          = ch_full_aln
    core_aln          = ch_core_aln
    core_vci          = ch_core_vci
}
