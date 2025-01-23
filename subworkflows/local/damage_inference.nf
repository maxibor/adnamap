include { PYDAMAGE_ANALYZE   } from '../../modules/nf-core/pydamage/analyze/main'
include { MAPDAMAGE2         } from '../../modules/nf-core/mapdamage2/main'

workflow DAMAGE_INFERENCE {
    take:
        input //meta, bam, bai, fasta, fai

    main:
        MAPDAMAGE2(
            input.map {
                meta, bam, bai, fasta, fai -> [meta, bam, fasta]
            }
        )

        PYDAMAGE_ANALYZE(
            input.map {
                meta, bam, bai, fasta, fai -> [meta, bam, bai]
            }
        )
}
