/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = Channel.fromPath(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.genomes) { ch_genomes = Channel.fromPath(params.genomes) } else { exit 1, 'Genomes sheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CREATE GENOMES CHANNEL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
ch_genomes
    .splitCsv(header:true, sep:',')
    .map { row -> [["genome_name": row.genome_name, "taxid": row.taxid], file(row.genome_path)] }
    .set { genomes }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CREATE DB CHANNEL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
Channel
    .fromPath(params.sam2lca_db)
    .first()
    .set { sam2lca_db }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK               } from '../subworkflows/local/input_check'
include { ALIGN_BOWTIE2             } from '../subworkflows/nf-core/align_bowtie2/main'
include { MERGE_SORT_INDEX_SAMTOOLS } from '../subworkflows/local/merge_sort_index_samtools'
include { VARIANT_CALLING           } from '../subworkflows/local/variant_calling'
include { BAM_SORT_SAMTOOLS         } from '../subworkflows/nf-core/bam_sort_samtools/main'
include { BAM_PICARD_MARKDUPLICATES } from '../subworkflows/local/bam_picard_markduplicates.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { SAMTOOLS_FAIDX                                   } from '../modules/nf-core/modules/samtools/faidx/main'
include { FASTQC as FASTQC_BEFORE ; FASTQC as FASTQC_AFTER } from '../modules/nf-core/modules/fastqc/main'
include { FASTP                                            } from '../modules/nf-core/modules/fastp/main'
include { GUNZIP                                           } from '../modules/nf-core/modules/gunzip/main'
include { BOWTIE2_BUILD                                    } from '../modules/nf-core/modules/bowtie2/build/main'
include { SAM2LCA                                          } from '../modules/local/sam2lca/main'
include { SAMTOOLS_INDEX as INDEX_PER_GENOME               } from '../modules/nf-core/modules/samtools/index/main'
include { BEDTOOLS_BAMTOBED                                } from '../modules/nf-core/modules/bedtools/bamtobed/main'
include { PRESEQ_LCEXTRAP                                  } from '../modules/nf-core/modules/preseq/lcextrap/main'
include { QUALIMAP_BAMQC                                   } from '../modules/nf-core/modules/qualimap/bamqc/main'
include { DAMAGEPROFILER                                   } from '../modules/nf-core/modules/damageprofiler/main'
include { MULTIQC                                          } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                      } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow ADNAMAP {

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FASTQ pre-processing
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC_BEFORE (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC_BEFORE.out.versions.first())

    FASTP (
        INPUT_CHECK.out.reads, params.save_trimmed_fail, params.save_merged
    )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    FASTQC_AFTER (
        FASTP.out.reads_merged.mix(FASTP.out.reads)
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FASTA pre-processing
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    genomes
        .branch {
            decompressed: it[1].toString().tokenize(".")[-1] != 'gz'
            compressed: it[1].toString().tokenize(".")[-1] == 'gz'
        }
        .set { genomes_fork}

    GUNZIP (
        genomes_fork.compressed
    )

    GUNZIP.out.gunzip
        .mix(genomes_fork.decompressed)
        .set { genomes_pre_processed }

    SAMTOOLS_FAIDX (
        genomes_pre_processed
    )

    BOWTIE2_BUILD (
        genomes_pre_processed
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Mixing read and reference channel
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */


    FASTP.out.reads_merged.mix(FASTP.out.reads) // meta_reads, merged_reads
        .combine(BOWTIE2_BUILD.out.index) // meta_genome, genome_index
        .map {
            meta_reads, reads, meta_genome, genome_index ->
                [
                    ['id': meta_reads.id + "_" + meta_genome.genome_name,
                     'genome_name': meta_genome.genome_name,
                     'taxid': meta_genome.taxid,
                     'sample_name': meta_reads.id,
                     'single_end': meta_reads.single_end

                    ],
                    reads,
                    genome_index
                ]
                // [meta_reads + meta_genome , reads, genome_index]
        }
        .set { ch_reads_genomes }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Alignment
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ALIGN_BOWTIE2 (
        ch_reads_genomes
    )
    ch_versions = ch_versions.mix(ALIGN_BOWTIE2.out.versions.first())

    ALIGN_BOWTIE2.out.bam.join(
        ALIGN_BOWTIE2.out.bai
    ).map {
        meta, bam, bai -> [['id':meta.sample_name], bam] // meta.id, bam
    }.groupTuple()
    .set { bams_synced }

    MERGE_SORT_INDEX_SAMTOOLS (
        bams_synced
    )

    SAM2LCA (
        MERGE_SORT_INDEX_SAMTOOLS.out.bam.join(
            MERGE_SORT_INDEX_SAMTOOLS.out.bai
        ),
        sam2lca_db
    )

    SAM2LCA.out.bam
    .transpose()
    .map {
        meta, bam ->
            def new_meta = [:]
                new_meta['id'] = meta.id
                new_meta['taxid'] = bam.baseName.toString().split("_taxid_")[-1].tokenize(".")[0]
            [new_meta['id'], new_meta['taxid'] , bam]
    }
    .join(
        ch_reads_genomes.map {
            meta, reads, genome_index -> [meta.sample_name, meta.taxid, meta]
        }, by: [0, 1]
    ).map {
        sample_name, taxid, bam, meta -> [meta, bam]
    }
    .set{ bam_split_by_ref} // meta, bam

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Duplicate Removal and Separation
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    BAM_SORT_SAMTOOLS {
        bam_split_by_ref
    }

    if (params.deduplicate == 'mark') {
        BAM_PICARD_MARKDUPLICATES (
            BAM_SORT_SAMTOOLS.out.bam
        )
        ch_versions = ch_versions.mix(BAM_PICARD_MARKDUPLICATES.out.versions.first())

        BEDTOOLS_BAMTOBED {
            BAM_PICARD_MARKDUPLICATES.out.bam
        }
        ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions.first())

        ch_bed_for_preseq = BEDTOOLS_BAMTOBED.out.bed.map{
            meta, bed ->
                def new_meta = meta.clone()
                new_meta['single_end'] = true

            [new_meta, bed]
        }

        PRESEQ_LCEXTRAP {
            ch_bed_for_preseq.dump()
        }
        ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions.first())

        sorted_bam_ch = BAM_PICARD_MARKDUPLICATES.out.bam

    } else {
        sorted_bam_ch = BAM_SORT_SAMTOOLS.out.bam
    }



    sorted_bam_ch.join(
        BAM_SORT_SAMTOOLS.out.bai
    ).map {
        meta, bam, bai -> [meta.taxid, meta, bam, bai] // taxid, meta, bam, bai
    }.combine(
        genomes_pre_processed
            .map{
                meta, fasta -> [meta.taxid, fasta]//taxid, fasta
            }
    , by: 0).combine(                      // taxid, meta, bam, bai, fasta
        SAMTOOLS_FAIDX.out.fai
            .map{
                meta, fai -> [meta.taxid, fai] // taxid, fai
            }
    , by: 0).map{ //taxid, meta, bam, bai, fasta, fai
        taxid, meta, bam, bai, fasta, fai -> [meta, bam, bai, fasta, fai] // meta, bam, bai, fasta, fai
    }.set {
        synced_ch
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Profiling and Genotyping
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */


    DAMAGEPROFILER (
        synced_ch
    )

    ch_versions = ch_versions.mix(DAMAGEPROFILER.out.versions.first())

    QUALIMAP_BAMQC (
        synced_ch
            .map{ meta, bam, bai, fasta, fai  -> [meta, bam] } // meta, bam
    )
    ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions.first())

    VARIANT_CALLING (
        synced_ch
    )

    ch_versions = ch_versions.mix(VARIANT_CALLING.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Post-processing
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowAdnamap.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_BEFORE.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_AFTER.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN_BOWTIE2.out.log_out.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_BAMQC.out.results.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(DAMAGEPROFILER.out.results.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(VARIANT_CALLING.out.stats.collect{it[1]}.ifEmpty([]))
    if (params.deduplicate == 'mark') {
        ch_multiqc_files = ch_multiqc_files.mix(BAM_PICARD_MARKDUPLICATES.out.metrics.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PRESEQ_LCEXTRAP.out.lc_extrap.collect{it[1]}.ifEmpty([]))
    }

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
