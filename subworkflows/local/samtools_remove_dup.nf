include { SAMTOOLS_VIEW                       } from  '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_FIXMATE                    } from '../../modules/nf-core/samtools/fixmate/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_COORD ;
          SAMTOOLS_SORT as SAMTOOLS_SORT_NAME } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_MARKDUP                    } from  '../../modules/nf-core/samtools/markdup/main'
include { SAMTOOLS_INDEX                      } from '../../modules/nf-core/samtools/index/main'

workflow SAMTOOLS_REMOVE_DUP {
    take:
        input //meta, bam, bai

    main:
        SAMTOOLS_SORT_NAME {
            input.map {
                    meta, bam, bai -> [meta, bam]
                }
        }

        SAMTOOLS_FIXMATE(
            SAMTOOLS_SORT_NAME.out.bam
        )

        SAMTOOLS_SORT_COORD(
            SAMTOOLS_FIXMATE.out.bam
        )

        SAMTOOLS_MARKDUP(
                SAMTOOLS_SORT_COORD.out.bam,
                []
            )

        SAMTOOLS_INDEX (
            SAMTOOLS_MARKDUP.out.bam
        )

        SAMTOOLS_VIEW(
            SAMTOOLS_MARKDUP.out.bam.join (
                SAMTOOLS_INDEX.out.bai
            ),
            [[],[]],
            []
        )

    emit:
        bam = SAMTOOLS_VIEW.out.bam
}
