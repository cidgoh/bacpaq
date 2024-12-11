include { SNIPPY_RUN        } from '../../modules/nf-core/snippy/run'
include { SNIPPY_CORE       } from '../../modules/nf-core/snippy/core'
include { MEDAKA_VARIANT    } from '../../modules/local/medaka/variant'
include { NUCMER            } from '../../modules/nf-core/nucmer'
include { DELTA2VCF         } from '../../modules/local/delta2vcf'


workflow VARIANT_CALLING {

    take:
        ch_contigs
        ch_reference_genome
        ch_long_reads

    main:

    ch_versions = Channel.empty()

    ch_reference_genome
        .combine(ch_contigs)
        .map {
            ref, contigs ->
            def meta        = [:]
                meta.id     = "testing"
                meta.ref    = "ref"
            [meta, ref, contigs] }
        .set{input_channel}

    NUCMER(input_channel)
    ch_versions = ch_versions.mix(NUCMER.out.versions)

    DELTA2VCF(NUCMER.out.delta)
    ch_versions = ch_versions.mix(DELTA2VCF.out.versions)

    ch_long_reads
        .map {
            fastq ->
            def meta        = [:]
                meta.id     = "testing"
                meta.ref    = "ref"
            [meta, fastq] }
        .set{input_long}

    MEDAKA_VARIANT(input_long, ch_reference_genome.first())
    MEDAKA_VARIANT.out.alignment_vcf.view()


    emit:
    versions        = ch_versions
    vcf_contigs     = DELTA2VCF.out.vcf
}
