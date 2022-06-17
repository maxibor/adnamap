include { SAMTOOLS_MERGE     } from  '../../../modules/nf-core/modules/samtools/merge/main'
include { BAM_SORT_SAMTOOLS  } from '../bam_sort_samtools/main'

workflow MERGE_SORT_INDEX_SAMTOOLS {
    take:
        ch_bam // [val(meta), [bams]]
    main:

    ch_versions = Channel.empty()

    SAMTOOLS_MERGE(ch_bam)
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    BAM_SORT_SAMTOOLS ( SAMTOOLS_MERGE.out.bam )
    ch_versions = ch_versions.mix(BAM_SORT_SAMTOOLS.out.versions)

    emit:
    bam = BAM_SORT_SAMTOOLS.out.bam
    bai = BAM_SORT_SAMTOOLS.out.bai

    versions = ch_versions
}
