include { IQTREE                                                     } from '../../modules/nf-core/iqtree'
include { IGVREPORTS as IGVREPORTS_VCF; IGVREPORTS as IGVREPORTS_BAM } from '../../modules/nf-core/igvreports'
include { SAMTOOLS_FAIDX                                             } from '../../modules/nf-core/samtools/faidx'
include { FAI2BED                                                    } from '../../modules/local/fai2bed/main'


workflow VARIANT_VIS {

    take:
        ch_vcf // [meta, vcf]
        ch_bam // [meta, bam]
        ch_bam_bai // [meta, bam_bai]
        ch_aln_fa // [meta, aln_fa]

    main:
    // initialize channels
    ch_versions        = Channel.empty()
    ch_igvreports_bam  = Channel.empty()
    ch_igvreports_vcf  = Channel.empty()
    ch_iqtree_tree     = Channel.empty()
    ch_iqtree_report   = Channel.empty()

    if (!params.skip_igvreports) {
        // index reference
        Channel
            .fromPath(params.reference_genome)
            .map { [ [id: 'reference'], it ]}
            .set { ch_reference_genome }
        SAMTOOLS_FAIDX(ch_reference_genome)
        ch_versions      = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        ch_reference_fai = SAMTOOLS_FAIDX.out.fai
        FAI2BED(ch_reference_fai)
        ch_reference_bed = FAI2BED.out.bed.map { meta, bed -> 
            [[ id: 'igv' ], bed]
        }

        // IGVREPORTS
        ch_IGVREPORTS_ref = ch_reference_genome.combine(ch_reference_fai, by:0)
        ch_IGVREPORTS_vcf = ch_vcf
            .map { it[1] }
            .collect()
            .map { [[ id: 'igv' ], it ] }
            .combine(ch_reference_bed, by: 0)
            .map { meta, vcf, bed -> [ meta, bed, vcf, [] ] }
        ch_IGVREPORTS_bam = ch_bam
            .map { it[1] }
            .collect()
            .map { [[ id: 'igv' ], it ] }
        ch_IGVREPORTS_bam_bai = ch_bam_bai
            .map { it[1] }
            .collect()
            .map { [[ id: 'igv' ], it ] }
        ch_IGVREPORTS_BAM_DATA=ch_IGVREPORTS_bam
            .combine(ch_IGVREPORTS_bam_bai, by: 0)
            .combine(ch_reference_bed, by: 0)
            .map { meta, bam, bai, bed -> [ meta, bed, bam, bai ] }

        IGVREPORTS_VCF(ch_IGVREPORTS_vcf, ch_IGVREPORTS_ref)
        ch_igvreports_vcf = IGVREPORTS_VCF.out.report
        ch_versions      = ch_versions.mix(IGVREPORTS_VCF.out.versions)
        IGVREPORTS_BAM(ch_IGVREPORTS_BAM_DATA, ch_IGVREPORTS_ref)
        ch_igvreports_bam = IGVREPORTS_BAM.out.report
        ch_versions      = ch_versions.mix(IGVREPORTS_BAM.out.versions)
    }
    
    
    // IQTREE
    if (!params.skip_iqtree) {

        ch_core_snp_aln = ch_aln_fa
            .map { meta, aln_fa ->
                [meta, aln_fa, []] 
            }

        IQTREE(
            ch_core_snp_aln, // alignment
            [], // tree_te
            [], // lmclust
            [], // mdef
            [], // partitions_equal
            [], // partitions_proportional
            [], // partitions_unlinked
            [], // guide_tree
            [], // sitefreq_in
            [], // constraint_tree
            [], // trees_z
            [], // suptree
            []  // trees_rf
        )
        ch_iqtree_tree     = IQTREE.out.phylogeny
        ch_iqtree_report   = IQTREE.out.report
        ch_versions        = ch_versions.mix(IQTREE.out.versions)
    }    

    emit:
    versions        = ch_versions
    igvreports_vcf  = ch_igvreports_vcf
    igvreports_bam  = ch_igvreports_bam
    iqtree_nwk      = ch_iqtree_tree
    iqtree_report   = ch_iqtree_report
}
