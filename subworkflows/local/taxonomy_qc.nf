//
// Importing the modules required for the sub-workflow
//


include { KRAKEN2_KRAKEN2           } from '../../modules/nf-core/kraken2/kraken2/main'
include { KRAKENTOOLS_KREPORT2KRONA } from '../../modules/nf-core/krakentools/kreport2krona/main'
include { BRACKEN_BRACKEN           } from '../../modules/nf-core/bracken/bracken/main'
include { MINIMAP2_ALIGN            } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_INDEX            } from '../../modules/nf-core/minimap2/index/main'
include { BWA_MEM                   } from '../../modules/nf-core/bwa/mem/main'
include { BWA_INDEX                 } from '../../modules/nf-core/bwa/index/main'
include { CENTRIFUGE_CENTRIFUGE     } from '../../modules/nf-core/centrifuge/centrifuge/main'
include { CENTRIFUGE_KREPORT        } from '../../modules/nf-core/centrifuge/kreport/main'
include { KRONA_KTUPDATETAXONOMY    } from '../../modules/nf-core/krona/ktupdatetaxonomy/main'
include { KRONA_KTIMPORTTAXONOMY    } from '../../modules/nf-core/krona/ktimporttaxonomy/main'
include { KRONA_KTIMPORTTEXT        } from '../../modules/nf-core/krona/ktimporttext/main'
include { SAMTOOLS_INDEX            } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FLAGSTAT         } from '../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_VIEW             } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_FASTQ            } from '../../modules/nf-core/samtools/fastq/main'

workflow TAXONOMY_QC {
    take:
    ch_reads_taxonomy
    reference_genome

    main:
    ch_versions = Channel.empty()
    ch_dehosted_reads = Channel.empty()
    kraken_report = Channel.empty()
    bracken_report = Channel.empty()

    if (params.classifier == "centrifuge") {
        ch_centrifuge_db = Channel.from(params.centrifuge_db)
        CENTRIFUGE_CENTRIFUGE(
            ch_reads_taxonomy,
            ch_centrifuge_db,
            params.save_unaligned,
            params.save_aligned
        )
        ch_versions = ch_versions.mix(CENTRIFUGE_CENTRIFUGE.out.versions)
        ch_tax_reads = CENTRIFUGE_CENTRIFUGE.out.fastq_mapped
        ch_tax_qc_unaligned_reads = CENTRIFUGE_CENTRIFUGE.out.fastq_unmapped
        CENTRIFUGE_KREPORT(
            CENTRIFUGE_CENTRIFUGE.out.report,
            ch_centrifuge_db
        )
        KRONA_KTUPDATETAXONOMY()
        ch_versions = ch_versions.mix(CENTRIFUGE_KREPORT.out.versions)
        KRONA_KTIMPORTTAXONOMY(
            CENTRIFUGE_KREPORT.out.kreport,
            KRONA_KTUPDATETAXONOMY.out.db
        )
        ch_versions = ch_versions.mix(KRONA_KTIMPORTTAXONOMY.out.versions)
    }
    else {
        unclassified_reads = params.unclassified_reads
        classified_reads = params.classified_reads
        //ch_kraken2_db=Channel.from(params.kraken2_db)
        if (!params.skip_kraken2) {
            KRAKEN2_KRAKEN2(
                ch_reads_taxonomy,
                params.kraken2_db,
                classified_reads,
                unclassified_reads
            )
            ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions)
            ch_tax_reads = KRAKEN2_KRAKEN2.out.classified_reads_fastq
            ch_tax_qc_unaligned_reads = KRAKEN2_KRAKEN2.out.unclassified_reads_fastq
            kraken_report = KRAKEN2_KRAKEN2.out.report
            if (!params.skip_kreport2krona) {
                KRAKENTOOLS_KREPORT2KRONA(
                    KRAKEN2_KRAKEN2.out.report
                )
                ch_versions = ch_versions.mix(KRAKENTOOLS_KREPORT2KRONA.out.versions)
                KRONA_KTUPDATETAXONOMY()
                ch_versions = ch_versions.mix(KRONA_KTUPDATETAXONOMY.out.versions)
                KRONA_KTIMPORTTEXT(
                    KRAKENTOOLS_KREPORT2KRONA.out.txt
                )
                ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT.out.versions)
            }
            if (!params.skip_bracken) {
                BRACKEN_BRACKEN(
                    KRAKEN2_KRAKEN2.out.report,
                    params.kraken2_db
                )
                bracken_report = BRACKEN_BRACKEN.out.reports
                ch_versions = ch_versions.mix(BRACKEN_BRACKEN.out.versions)
            }
        }
    }
    ch_tax_qc_reads = ch_reads_taxonomy
    if (!params.skip_dehosting) {
        if (params.dehosting_aligner == 'bwa') {

            BWA_INDEX(
                [[id: params.ref_genome_id], reference_genome]
            )
            ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
            ref_genome = [[id: params.ref_genome_id], reference_genome]
            reads = ch_reads_taxonomy
            sorted = params.bwa_sort_bam
            BWA_MEM(
                reads,
                BWA_INDEX.out.index,
                ref_genome,
                sorted
            )
            ch_versions = ch_versions.mix(BWA_MEM.out.versions)
            ch_mapped_bam = BWA_MEM.out.bam
        }
        else {
            ref_genome = [[id: params.ref_genome_id], reference_genome]
            reads = ch_reads_taxonomy
            // bam_format = params.bam_format
            // true
            // cigar_paf_format = params.cigar_paf_format
            // false
            // cigar_bam = params.cigar_bam
            // false
            MINIMAP2_ALIGN(
                reads,
                ref_genome,
                params.bam_format,
                'bai',
                params.cigar_paf_format,
                params.cigar_bam
            )
            ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
            ch_mapped_bam = MINIMAP2_ALIGN.out.bam
        }

        SAMTOOLS_INDEX(
            ch_mapped_bam
        )
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
        bam_index = SAMTOOLS_INDEX.out.bai

        SAMTOOLS_FLAGSTAT(
            ch_mapped_bam.combine(bam_index, by: 0)
        )
        ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

        SAMTOOLS_VIEW(
            ch_mapped_bam.combine(bam_index, by: 0),
            [[], []],
            []
        )
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)

        SAMTOOLS_FASTQ(
            SAMTOOLS_VIEW.out.bam,
            params.interleaved
        )
        ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions)

        // illumina reads
        if (params.interleaved) {
            ch_dehosted_reads = ch_dehosted_reads.mix(SAMTOOLS_FASTQ.out.interleaved.filter { it[0].mode == 'illumina' })
        }
        else {
            ch_dehosted_reads = ch_dehosted_reads.mix(SAMTOOLS_FASTQ.out.fastq.filter { it[0].mode == 'illumina' })
        }
        // nanopore reads
        ch_dehosted_reads = ch_dehosted_reads.mix(SAMTOOLS_FASTQ.out.other.filter { it[0].mode == 'nanopore' })
        // output dehosted reads
        ch_tax_qc_reads = ch_dehosted_reads
    }

    emit:
    reads          = ch_tax_qc_reads // channel: [ val(meta), [ reads ] ]
    kraken_report
    bracken_report
    versions       = ch_versions // channel: [ versions.yml ]
}
