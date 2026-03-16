include { SNIPPY_RUN   as SNIPPY_RUN_READS  } from '../../modules/nf-core/snippy/run'
include { SNIPPY_RUN_CONTIGS } from '../../modules/local/snippy/run'
include { SNIPPY_CORE    } from '../../modules/nf-core/snippy/core'
include { GUBBINS        } from '../../modules/nf-core/gubbins'
include { MEDAKA_VARIANT } from '../../modules/local/medaka/variant'
include { DELTA2VCF      } from '../../modules/local/delta2vcf'
include { CORESNPFILTER  } from '../../modules/local/coresnpfilter'
include { SAMTOOLS_VIEW  } from '../../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index'
include { SAMTOOLS_SORT  } from '../../modules/nf-core/samtools/sort'
include { TABIX_TABIX as TABIX_SNIPPY_READS; TABIX_TABIX as TABIX_MEDAKA; TABIX_TABIX as TABIX_SNIPPY_CONTIGS; TABIX_TABIX as TABIX_CORE } from '../../modules/nf-core/tabix/tabix'
include { TABIX_BGZIP as BGZIP_SNIPPY_READS; TABIX_BGZIP as BGZIP_MEDAKA; TABIX_BGZIP as BGZIP_SNIPPY_CONTIGS; TABIX_BGZIP as BGZIP_CORE } from '../../modules/nf-core/tabix/bgzip'


workflow VARIANT_CALLING {
    take:
    ch_reads  // channel containing path to sequence data in FASTQ format
    ch_genome // channel containing path to reference genome in FASTA format

    main:
    // initialize channels
    ch_versions = Channel.empty()
    ch_vcf_snippy_contigs = Channel.empty()
    ch_vcf_snippy_reads = Channel.empty()
    ch_vcf_medaka = Channel.empty()
    ch_vcf_gubbins = Channel.empty()
    ch_vcf_bgz_snippy_contigs = Channel.empty()
    ch_vcf_bgz_snippy_reads = Channel.empty()
    ch_vcf_bgz_medaka = Channel.empty()
    ch_vcf_gubbins = Channel.empty()
    ch_sam_snippy_contigs = Channel.empty()
    ch_bam_snippy_contigs = Channel.empty()
    ch_bam_snippy_reads = Channel.empty()
    ch_bam_medaka = Channel.empty()
    ch_bai_snippy_reads = Channel.empty()
    ch_bai_medaka = Channel.empty()
    ch_bai_snippy_contigs = Channel.empty()
    ch_vci_snippy_reads = Channel.empty()
    ch_vci_medaka = Channel.empty()
    ch_vci_snippy_contigs = Channel.empty()
    // snippy output channels
    ch_tab_snippy_reads = Channel.empty()
    ch_csv_snippy_reads = Channel.empty()
    ch_html_snippy_reads = Channel.empty()
    ch_bed_snippy_reads = Channel.empty()
    ch_gff_snippy_reads = Channel.empty()
    ch_log_snippy_reads = Channel.empty()
    ch_aligned_fa_snippy_reads = Channel.empty()
    ch_consensus_fa_snippy_reads = Channel.empty()
    ch_consensus_subs_fa_snippy_reads = Channel.empty()
    ch_raw_vcf_snippy_reads = Channel.empty()
    ch_filt_vcf_snippy_reads = Channel.empty()
    ch_vcf_csi_snippy_reads = Channel.empty()
    ch_vcf_gz_snippy_reads = Channel.empty()
    ch_vcf_csi_snippy_reads = Channel.empty()
    ch_vcf_gz_snippy_reads = Channel.empty()
    ch_txt_snippy_reads = Channel.empty()
    ch_txt_snippy_contigs = Channel.empty()
    
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
    ch_genome.map{
        meta, path -> [meta, path[0]]
        }
        .set { ch_snippy_contig }

    // GENOME: RUN SNIPPY CONTIG
    SNIPPY_RUN_CONTIGS(ch_snippy_contig, ch_reference_genome.first())

    ch_vcf_snippy_contigs = SNIPPY_RUN_CONTIGS.out.vcf
    ch_bam_snippy_contigs = SNIPPY_RUN_CONTIGS.out.bam
    ch_bai_snippy_contigs = SNIPPY_RUN_CONTIGS.out.bai
    ch_tab_snippy_contigs = SNIPPY_RUN_CONTIGS.out.tab
    ch_csv_snippy_contigs = SNIPPY_RUN_CONTIGS.out.csv
    ch_html_snippy_contigs = SNIPPY_RUN_CONTIGS.out.html
    ch_bed_snippy_contigs = SNIPPY_RUN_CONTIGS.out.bed
    ch_gff_snippy_contigs = SNIPPY_RUN_CONTIGS.out.gff
    ch_log_snippy_contigs = SNIPPY_RUN_CONTIGS.out.log
    ch_aligned_fa_snippy_contigs = SNIPPY_RUN_CONTIGS.out.aligned_fa
    ch_consensus_fa_snippy_contigs = SNIPPY_RUN_CONTIGS.out.consensus_fa
    ch_consensus_subs_fa_snippy_contigs = SNIPPY_RUN_CONTIGS.out.consensus_subs_fa
    ch_raw_vcf_snippy_contigs = SNIPPY_RUN_CONTIGS.out.raw_vcf
    ch_filt_vcf_snippy_contigs = SNIPPY_RUN_CONTIGS.out.filt_vcf
    ch_vcf_csi_snippy_contigs = SNIPPY_RUN_CONTIGS.out.vcf_csi
    ch_vcf_gz_snippy_contigs = SNIPPY_RUN_CONTIGS.out.vcf_gz
    ch_vcf_csi_snippy_contigs = SNIPPY_RUN_CONTIGS.out.vcf_csi
    ch_vcf_gz_snippy_contigs = SNIPPY_RUN_CONTIGS.out.vcf_gz
    ch_txt_snippy_contigs = SNIPPY_RUN_CONTIGS.out.txt
    
    BGZIP_SNIPPY_CONTIGS(ch_vcf_snippy_contigs) // compress the VCF file
    ch_versions = ch_versions.mix(BGZIP_SNIPPY_CONTIGS.out.versions)
    ch_vcf_bgz_snippy_contigs = BGZIP_SNIPPY_CONTIGS.out.output
    TABIX_SNIPPY_CONTIGS(ch_vcf_bgz_snippy_contigs) // index the VCF file
    ch_versions = ch_versions.mix(TABIX_SNIPPY_CONTIGS.out.versions)
    ch_vci_snippy_contigs = TABIX_SNIPPY_CONTIGS.out.tbi


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
    ch_snippy_reads = ch_input_seq.illumina
    SNIPPY_RUN_READS(ch_snippy_reads, ch_reference_genome.first())
    ch_vcf_snippy_reads = SNIPPY_RUN_READS.out.vcf
    ch_bam_snippy_reads = SNIPPY_RUN_READS.out.bam
    ch_bai_snippy_reads = SNIPPY_RUN_READS.out.bai
    ch_tab_snippy_reads = SNIPPY_RUN_READS.out.tab
    ch_csv_snippy_reads = SNIPPY_RUN_READS.out.csv
    ch_html_snippy_reads = SNIPPY_RUN_READS.out.html
    ch_bed_snippy_reads = SNIPPY_RUN_READS.out.bed
    ch_gff_snippy_reads = SNIPPY_RUN_READS.out.gff
    ch_log_snippy_reads = SNIPPY_RUN_READS.out.log
    ch_aligned_fa_snippy_reads = SNIPPY_RUN_READS.out.aligned_fa
    ch_consensus_fa_snippy_reads = SNIPPY_RUN_READS.out.consensus_fa
    ch_consensus_subs_fa_snippy_reads = SNIPPY_RUN_READS.out.consensus_subs_fa
    ch_raw_vcf_snippy_reads = SNIPPY_RUN_READS.out.raw_vcf
    ch_filt_vcf_snippy_reads = SNIPPY_RUN_READS.out.filt_vcf
    ch_vcf_csi_snippy_reads = SNIPPY_RUN_READS.out.vcf_csi
    ch_vcf_gz_snippy_reads = SNIPPY_RUN_READS.out.vcf_gz
    ch_vcf_csi_snippy_reads = SNIPPY_RUN_READS.out.vcf_csi
    ch_vcf_gz_snippy_reads = SNIPPY_RUN_READS.out.vcf_gz
    ch_txt_snippy_reads = SNIPPY_RUN_READS.out.txt
    

    BGZIP_SNIPPY_READS(ch_vcf_snippy_reads) // compress the VCF file
    ch_versions = ch_versions.mix(BGZIP_SNIPPY_READS.out.versions)
    ch_vcf_bgz_snippy_reads = BGZIP_SNIPPY_READS.out.output
    TABIX_SNIPPY_READS(ch_vcf_bgz_snippy_reads) // index the VCF file
    ch_versions = ch_versions.mix(TABIX_SNIPPY_READS.out.versions)
    ch_vci_snippy_reads = TABIX_SNIPPY_READS.out.tbi

    ch_snippy_reads_aligned_fa = SNIPPY_RUN_READS.out.aligned_fa
        .map { it[1] }
        .collect()
        .map { [[id: 'core_aln'], it] }

    ch_snippy_reads_vcf = SNIPPY_RUN_READS.out.vcf
        .map { it[1] }
        .collect()
        .map { [[id: 'core_aln'], it] }

    // generate core SNP alignment
    ch_snippy_reads_core = ch_snippy_reads_vcf.combine(ch_snippy_reads_aligned_fa, by: 0)
    SNIPPY_CORE(ch_snippy_reads_core, ch_reference_genome)
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
    vcf_snippy_contigs        = ch_vcf_snippy_contigs
    vcf_snippy_reads        = ch_vcf_snippy_reads
    vcf_medaka        = ch_vcf_medaka
    vcf_gubbins       = ch_vcf_gubbins
    vcf_bgz_snippy_contigs        = ch_vcf_bgz_snippy_contigs
    vcf_bgz_snippy_reads     = ch_vcf_bgz_snippy_reads
    vcf_bgz_medaka        = ch_vcf_bgz_medaka
    bam_snippy_reads        = ch_bam_snippy_reads
    bam_medaka        = ch_bam_medaka
    bam_snippy_contigs        = ch_bam_snippy_contigs
    bai_snippy_reads        = ch_bai_snippy_reads
    bai_medaka        = ch_bai_medaka
    bai_snippy_contigs        = ch_bai_snippy_contigs
    vci_snippy_reads        = ch_vci_snippy_reads
    vci_medaka        = ch_vci_medaka
    vci_snippy_contigs        = ch_vci_snippy_contigs

    // snippy reads outputs
    tab_snippy_reads        = ch_tab_snippy_reads
    csv_snippy_reads        = ch_csv_snippy_reads
    html_snippy_reads       = ch_html_snippy_reads
    bed_snippy_reads        = ch_bed_snippy_reads
    gff_snippy_reads        = ch_gff_snippy_reads
    log_snippy_reads        = ch_log_snippy_reads
    aligned_fa_snippy_reads = ch_aligned_fa_snippy_reads
    consensus_fa_snippy_reads = ch_consensus_fa_snippy_reads
    consensus_subs_fa_snippy_reads = ch_consensus_subs_fa_snippy_reads
    raw_vcf_snippy_reads    = ch_raw_vcf_snippy_reads
    filt_vcf_snippy_reads   = ch_filt_vcf_snippy_reads
    vcf_csi_snippy_reads    = ch_vcf_csi_snippy_reads
    vcf_gz_snippy_reads     = ch_vcf_gz_snippy_reads
    txt_snippy_reads     = ch_txt_snippy_reads

    // snippy contigs outputs
    tab_snippy_contigs        = ch_tab_snippy_contigs
    csv_snippy_contigs        = ch_csv_snippy_contigs
    html_snippy_contigs       = ch_html_snippy_contigs
    bed_snippy_contigs        = ch_bed_snippy_contigs
    gff_snippy_contigs        = ch_gff_snippy_contigs
    log_snippy_contigs        = ch_log_snippy_contigs
    aligned_fa_snippy_contigs = ch_aligned_fa_snippy_contigs
    consensus_fa_snippy_contigs = ch_consensus_fa_snippy_contigs
    consensus_subs_fa_snippy_contigs = ch_consensus_subs_fa_snippy_contigs
    raw_vcf_snippy_contigs    = ch_raw_vcf_snippy_contigs
    filt_vcf_snippy_contigs   = ch_filt_vcf_snippy_contigs
    vcf_csi_snippy_contigs    = ch_vcf_csi_snippy_contigs
    vcf_gz_snippy_contigs     = ch_vcf_gz_snippy_contigs
    txt_snippy_contigs     = ch_txt_snippy_contigs
    
    // snippy core outputs
    core_tab          = ch_core_tab
    core_vcf          = ch_core_vcf
    core_txt          = ch_core_txt
    full_aln          = ch_full_aln
    core_aln          = ch_core_aln
    core_vci          = ch_core_vci
}
