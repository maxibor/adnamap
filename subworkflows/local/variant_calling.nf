include { FREEBAYES          } from '../../modules/nf-core/modules/freebayes/main'
include { BCFTOOLS_INDEX     } from '../../modules/nf-core/modules/bcftools/index/main'
include { BCFTOOLS_VIEW      } from '../../modules/nf-core/modules/bcftools/view/main'
include { BCFTOOLS_STATS     } from '../../modules/nf-core/modules/bcftools/stats/main'
include { BCFTOOLS_CONSENSUS } from '../../modules/nf-core/modules/bcftools/consensus/main'


workflow VARIANT_CALLING {
    take:
    synced_ch // meta, bam, bai, fasta, fai

    main:
    ch_versions = Channel.empty()

    FREEBAYES (
        synced_ch
    )

    ch_versions = ch_versions.mix(FREEBAYES.out.versions)

    BCFTOOLS_INDEX (
        FREEBAYES.out.vcf
    )

    BCFTOOLS_STATS (
        FREEBAYES.out.vcf, []
    )

    BCFTOOLS_VIEW (
        FREEBAYES.out.vcf.join(BCFTOOLS_INDEX.out.tbi), [], [], []
    )

    BCFTOOLS_CONSENSUS (
        BCFTOOLS_VIEW.out.vcf
            .join(BCFTOOLS_INDEX.out.tbi)
            .join(synced_ch.map {it -> [it[0], it[3]]}) // meta, fasta
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions)

    emit:
    vcf= FREEBAYES.out.vcf
    tbi= BCFTOOLS_INDEX.out.tbi
    stats= BCFTOOLS_STATS.out.stats
    consensus= BCFTOOLS_CONSENSUS.out.fasta
    versions= ch_versions
}
