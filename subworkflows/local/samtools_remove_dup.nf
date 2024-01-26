include { SAMTOOLS_VIEW } from  '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_MARKDUP  } from  '../../modules/nf-core/samtools/markdup/main'

workflow SAMTOOLS_REMOVE_DUP {
    take:
        input //meta, bam, bai
    main:
        SAMTOOLS_MARKDUP(
                input.map {
                    meta, bam, bai -> [meta, bam]
                }
            )
        SAMTOOLS_VIEW(
            SAMTOOLS_MARKDUP.out.bam.join (
                input.map {
                    meta, bam, bai -> [meta, bai]
                }
            ),
            [],
            []
        )
    emit:
        bam = SAMTOOLS_VIEW.out.bam
}
