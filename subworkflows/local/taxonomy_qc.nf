//
// Importing the modules required for the sub-workflow
//


include { KRAKEN2_KRAKEN2               } from '../../modules/nf-core/kraken2/kraken2/main'
include { KRAKENTOOLS_KREPORT2KRONA     } from '../../modules/nf-core/krakentools/kreport2krona/main'
include { BRACKEN_BRACKEN               } from '../../modules/nf-core/bracken/bracken/main'
include { MINIMAP2_ALIGN as MINIMAP2_ILLUMINA               } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_NANOPORE               } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_INDEX                } from '../../modules/nf-core/minimap2/index/main'
include { BWA_MEM                       } from '../../modules/nf-core/bwa/mem/main'
include { BWA_INDEX                     } from '../../modules/nf-core/bwa/index/main'
include { CENTRIFUGE_CENTRIFUGE         } from '../../modules/nf-core/centrifuge/centrifuge/main'
include { CENTRIFUGE_KREPORT            } from '../../modules/nf-core/centrifuge/kreport/main'
include { KRONA_KTUPDATETAXONOMY        } from '../../modules/nf-core/krona/ktupdatetaxonomy/main'
include { KRONA_KTIMPORTTAXONOMY        } from '../../modules/nf-core/krona/ktimporttaxonomy/main'
include { KRONA_KTIMPORTTEXT            } from '../../modules/nf-core/krona/ktimporttext/main'
include { SAMTOOLS_INDEX                } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FLAGSTAT             } from '../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_VIEW                 } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_FASTQ                } from '../../modules/nf-core/samtools/fastq/main'

workflow TAXONOMY_QC {
    take:
    ch_reads_taxonomy 
    reference_genome

    main:
    ch_versions = Channel.empty()
    if (params.classifier=="centrifuge"){
        ch_centrifuge_db=Channel.from(params.centrifuge_db)
        CENTRIFUGE_CENTRIFUGE(
            ch_reads_taxonomy,
            ch_centrifuge_db,
            params.save_unaligned,
            params.save_aligned,
            params.sam_format
        )
        ch_tax_reads = CENTRIFUGE_CENTRIFUGE.out.fastq_mapped
        ch_tax_qc_unaligned_reads = CENTRIFUGE_CENTRIFUGE.out.fastq_unmapped
        CENTRIFUGE_KREPORT(
            CENTRIFUGE_CENTRIFUGE.out.report,
            ch_centrifuge_db
        )
        KRONA_KTIMPORTTAXONOMY(
                    CENTRIFUGE_KREPORT.out.kreport
        )
    }
    else{
        unclassified_reads=params.unclassified_reads
        classified_reads=params.classified_reads
        //ch_kraken2_db=Channel.from(params.kraken2_db)
        if (!params.skip_kraken2) {
            KRAKEN2_KRAKEN2(
                ch_reads_taxonomy,
                params.kraken2_db,
                classified_reads,
                unclassified_reads
            )
            ch_tax_reads = KRAKEN2_KRAKEN2.out.classified_reads_fastq
            ch_tax_qc_unaligned_reads = KRAKEN2_KRAKEN2.out.unclassified_reads_fastq
            if (!params.skip_kreport2krona) {
                KRAKENTOOLS_KREPORT2KRONA(
                    KRAKEN2_KRAKEN2.out.report
                )
                KRONA_KTUPDATETAXONOMY(

                )
                KRONA_KTIMPORTTEXT(
                    KRAKENTOOLS_KREPORT2KRONA.out.txt
                )
            }
            if (!params.skip_bracken) {
                BRACKEN_BRACKEN(
                    KRAKEN2_KRAKEN2.out.report,
                    params.kraken2_db
                )
            }

        }
         
    }

    if (!params.skip_dehosting){
        if (params.dehosting_aligner=='bwa') {

            BWA_INDEX(
                [[id: params.ref_genome_id], reference_genome]
            )
            ref_genome = BWA_INDEX.out.index
            reads = ch_reads_taxonomy
            sorted=params.bwa_sort_bam
            BWA_MEM(
                reads,
                ref_genome,
                sorted
            )
            ch_mapped_bam=BWA_MEM.out.bam
        }
        else {
            ref_genome = reference_genome
            reads = ch_reads_taxonomy
            bam_format  = params.bam_format //true
            cigar_paf_format = params.cigar_paf_format //false
            cigar_bam = params.cigar_bam //false
            if(params.mode == 'nanopore'){
                MINIMAP2_NANOPORE ( 
                reads, 
                ref_genome, 
                bam_format, 
                cigar_paf_format, 
                cigar_bam 
                )
            ch_mapped_bam=MINIMAP2_NANOPORE.out.bam
            }
            else{
                MINIMAP2_ILLUMINA ( 
                reads, 
                ref_genome, 
                bam_format, 
                cigar_paf_format, 
                cigar_bam 
                )
            ch_mapped_bam=MINIMAP2_ILLUMINA.out.bam
            } 
        }

        SAMTOOLS_INDEX(
            ch_mapped_bam
            )
        bam_index = SAMTOOLS_INDEX.out.bai

        SAMTOOLS_FLAGSTAT(
            ch_mapped_bam.combine(bam_index, by: 0)
            )

        SAMTOOLS_VIEW(
            ch_mapped_bam.combine(bam_index, by: 0),
            [],
            [])
        interleaved=params.interleaved
        SAMTOOLS_FASTQ(
             SAMTOOLS_VIEW.out.bam, interleaved)
        if (!params.mode == "nanopore"){
            if (interleaved == false){
                ch_dehosted_reads = SAMTOOLS_FASTQ.out.fastq
            }
            else{
                ch_dehosted_reads = SAMTOOLS_FASTQ.out.interleaved
            }
        }
        else{
            ch_dehosted_reads = SAMTOOLS_FASTQ.out.other
        }
        
    }
    
    emit:
    ch_tax_qc_reads = ch_dehosted_reads                              // channel: [ val(meta), [ reads ] ]
    //ch_tax_unaligned_reads = ch_tax_qc_unaligned_reads                    // channel: [ val(meta), [ reads ] ]
    //versions = TAXONOMY_QC.out.versions                                  // channel: [ versions.yml ]
}
