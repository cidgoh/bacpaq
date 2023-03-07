//
// Importing the modules required for the sub-workflow
//


include { KRAKEN2_KRAKEN2 } from '../../modules/nf-core/kraken2/kraken2/main'
include { KRAKENTOOLS_COMBINEKREPORTS } from '../../modules/nf-core/krakentools/combinekreports/main'
include { KRAKENTOOLS_KREPORT2KRONA } from '../../modules/nf-core/krakentools/kreport2krona/main'
include { BRACKEN_BRACKEN } from '../../modules/nf-core/bracken/bracken/main'
include { BRACKEN_COMBINEBRACKENOUTPUTS } from '../../modules/nf-core/bracken/combinebrackenoutputs/main'
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main'
include { BWA_MEM } from '../../modules/nf-core/bwa/mem/main'
include { BWA_INDEX } from '../../modules/nf-core/bwa/index/main'
include { CENTRIFUGE_CENTRIFUGE } from '../../modules/nf-core/centrifuge/centrifuge/main'
include { CENTRIFUGE_KREPORT } from '../../modules/nf-core/centrifuge/kreport/main'
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_FLAGSTAT } from '../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_VIEW } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_FASTQ } from '../../modules/nf-core/samtools/fastq/main'

workflow TAXONOMY_QC {
    take:
    ch_reads_taxonomy
    classified_reads
    unclassified_reads 
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
        CENTRIFUGE_KREPORT(
            CENTRIFUGE_CENTRIFUGE.out.report,
            ch_centrifuge_db
        )
    }
    else{
        ch_kraken2_db=Channel.from(params.kraken2_db)
        if (!params.skip_kraken2) {
            KRAKEN2_KRAKEN2(
                ch_reads_taxonomy,
                ch_kraken2_db,
                classified_reads,
                unclassified_reads
            )
            if (!params.skip_combinekreports) {
                KRAKENTOOLS_COMBINEKREPORTS(
                    KRAKEN2_KRAKEN2.out.report.collect()
                )
                    
            }
            if (!params.skip_kreport2krona) {
                KRAKENTOOLS_KREPORT2KRONA(
                    KRAKEN2_KRAKEN2.out.report
                )
        
            }
            if (!params.skip_bracken) {
                BRACKEN_BRACKEN(
                    KRAKEN2_KRAKEN2.out.report,
                    ch_kraken2_db
                )
            }
            if (!params.skip_combinebrackenoutputs) {
                BRACKEN_COMBINEBRACKENOUTPUTS(
                    BRACKEN_BRACKEN.out.reports.collect()
                )
            }

        }
         
    }
    


    if (!params.skip_dehosting){
        if (params.dehosting_aligner=='bwa') {
                
        }
        else {
            // Minimap2
            
            reads = ch_reads_taxonomy
            bam_format  = params.bam_format //true
            cigar_paf_format = params.cigar_paf_format //false
            cigar_bam = params.cigar_bam //false
            MINIMAP2_ALIGN ( reads, ref_genome, bam_format, cigar_paf_format, cigar_bam )
            ch_mapped_bam=MINIMAP2_ALIGN.out.bam
        }

        SAMTOOLS_FLAGSTAT(
            ch_mapped_bam)

        SAMTOOLS_VIEW(
            ch_mapped_bam)

        SAMTOOLS_FASTQ(
            SAMTOOLS_VIEW.out.bam)
        // samtools
    }
    versions = TAXONOMY_QC.out.versions.first()
    
    emit:
    classified_reads                           // channel: [ val(meta), [ reads ] ]
    unclassified_reads                         // channel: [ val(meta), [ reads ] ]
    versions                                  // channel: [ versions.yml ]
}
