include { FREEBAYES } from '../../modules/nf-core/freebayes/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_PRE;
          BCFTOOLS_INDEX as BCFTOOLS_INDEX_POST
        } from '../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_VIEW      } from '../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_STATS     } from '../../modules/nf-core/bcftools/stats/main'
include { BCFTOOLS_CONSENSUS } from '../../modules/nf-core/bcftools/consensus/main'


workflow VARIANT_CALLING {
    take:
    synced_ch // meta, bam, bai, fasta, fai

    main:
    ch_versions = Channel.empty()

    FREEBAYES (
        synced_ch
    )

    ch_versions = ch_versions.mix(FREEBAYES.out.versions)

    BCFTOOLS_INDEX_PRE (
        FREEBAYES.out.vcf
    )

    BCFTOOLS_STATS (
        FREEBAYES.out.vcf.join(BCFTOOLS_INDEX_PRE.out.tbi), [], [], [], [], []
    )

    BCFTOOLS_VIEW (
        FREEBAYES.out.vcf.join(BCFTOOLS_INDEX_PRE.out.tbi), [], [], []
    )

    BCFTOOLS_INDEX_POST(
        BCFTOOLS_VIEW.out.vcf
    )

    BCFTOOLS_CONSENSUS (
        BCFTOOLS_VIEW.out.vcf
            .join(BCFTOOLS_INDEX_POST.out.tbi)
            .join(synced_ch.map {it -> [it[0], it[3]]}) // meta, fasta
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions)

    emit:
    vcf= FREEBAYES.out.vcf
    stats= BCFTOOLS_STATS.out.stats
    consensus= BCFTOOLS_CONSENSUS.out.fasta
    versions= ch_versions
}
