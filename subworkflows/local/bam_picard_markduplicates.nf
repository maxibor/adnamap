//
// Sort, Deduplicate and Resort
//

include { PICARD_MARKDUPLICATES } from '../../modules/nf-core/modules/picard/markduplicates/main'
include { SAMTOOLS_SORT         } from '../../modules/nf-core/modules/samtools/sort/main'
include { SAMTOOLS_INDEX        } from '../../modules/nf-core/modules/samtools/index/main'


workflow BAM_PICARD_MARKDUPLICATES {

    take:
    ch_bam  // channel: [ val(meta), [ bam ] ]

    main:

    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()


    PICARD_MARKDUPLICATES {
        ch_bam
    }

    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())

    SAMTOOLS_SORT {
        PICARD_MARKDUPLICATES.out.bam
    }
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_INDEX {
        SAMTOOLS_SORT.out.bam
    }
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    unsorted_bam = PICARD_MARKDUPLICATES.out.bam      // channel: [ val(meta), [ bam ] ]
    unsorted_bai = PICARD_MARKDUPLICATES.out.bai      // OPTIONAL channel: [ val(meta), [ bai ] ]
    metric       = PICARD_MARKDUPLICATES.out.metrics  // channel: [ val(meta), [ metrics ] ]
    bam          = SAMTOOLS_SORT.out.bam              // channel: [ val(meta), [ bam ] ]
    bai          = SAMTOOLS_INDEX.out.bai             // channel: [ val(meta), [ bai ] ]
    csi          = SAMTOOLS_INDEX.out.csi             // channel: [ val(meta), [ csi ] ]


    versions = ch_versions                            // channel: [ versions.yml ]

}
