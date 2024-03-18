/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

// // Validate input parameters
// WorkflowMag.initialise(params, log)

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
    .map { row -> [ ["genome_name": row.genome_name, "taxid": row.taxid, "ploidy": row.ploidy], file(row.genome_path), file(row.genome_index, type: 'dir') ] }
    .set { genomes }

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
include { SAM2LCA_DB                } from '../subworkflows/local/sam2lca_db'
include { SAMTOOLS_REMOVE_DUP       } from '../subworkflows/local/samtools_remove_dup'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { SAMTOOLS_FAIDX                                   } from '../modules/nf-core/samtools/faidx/main'
include { FASTQC as FASTQC_BEFORE ; FASTQC as FASTQC_AFTER } from '../modules/nf-core/fastqc/main'
include { FASTP                                            } from '../modules/nf-core/fastp/main'
include { GUNZIP as GUNZIP4IDX ; GUNZIP as GUNZIP4GENOME   } from '../modules/nf-core/gunzip/main'
include { UNTAR                                            } from '../modules/nf-core/untar/main'
include { BOWTIE2_BUILD                                    } from '../modules/nf-core/bowtie2/build/main'
include { PRESEQ_LCEXTRAP                                  } from '../modules/nf-core/preseq/lcextrap/main'
include { PRESEQ_CCURVE                                    } from '../modules/nf-core/preseq/ccurve/main'
include { PLOT_PRESEQ                                      } from '../modules/local/plot_preseq'
include { SAM2LCA                                          } from '../modules/local/sam2lca/main'
include { SAMTOOLS_INDEX as INDEX_PER_GENOME               } from '../modules/nf-core/samtools/index/main'
include { QUALIMAP_BAMQC                                   } from '../modules/nf-core/qualimap/bamqc/main'
include { MAPDAMAGE2                                       } from '../modules/nf-core/mapdamage2/main'
include { DAMAGEPROFILER                                   } from '../modules/nf-core/damageprofiler/main'
include { MULTIQC                                          } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                      } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { SAM2LCA_MERGE                                    } from '../modules/local/sam2lca_merge'
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
        INPUT_CHECK.out.reads, [], params.save_trimmed_fail, params.save_merged
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
            decompressed: it[1].toString().tokenize(".")[-1] != 'gz' && it[2].isEmpty()
            compressed: it[1].toString().tokenize(".")[-1] == 'gz' && it[2].isEmpty()
            has_index_decompressed: ! it[2].isEmpty() && it[2].toString().tokenize(".")[-1] != 'gz'
            has_index_compressed: ! it[2].isEmpty() && it[2].toString().tokenize(".")[-1] == 'gz'
        }
        .set { genomes_idx_fork }

    genomes
        .branch {
            decompressed: it[1].toString().tokenize(".")[-1] != 'gz' && ! it[2].isEmpty()
            compressed: it[1].toString().tokenize(".")[-1] == 'gz' && ! it[2].isEmpty()
        }
        .set { genomes_fasta_fork }

    GUNZIP4IDX (
        genomes_idx_fork.compressed
        .map {
            genome_meta, fasta, index -> [genome_meta, fasta]
        }
    )

    UNTAR(
        genomes_idx_fork.has_index_compressed.map {
            meta, genome_fasta, genome_index ->
            [meta, genome_index]
        }
    )

    GUNZIP4IDX.out.gunzip
        .mix(genomes_idx_fork.decompressed)
        .set { genomes_pre_processed }

    BOWTIE2_BUILD (
        genomes_pre_processed
    )

    ch_indices = BOWTIE2_BUILD.out.index.mix(
            genomes_idx_fork.has_index_decompressed
            .map {
                meta_genome, genome_fasta, genome_index ->
                [meta_genome, genome_index]
            }
        ).mix (
                UNTAR.out.untar
        )

    GUNZIP4GENOME (
        genomes_fasta_fork.compressed
        .map {
            genome_meta, fasta, index -> [genome_meta, fasta]
        }
    )

    GUNZIP4GENOME.out.gunzip.mix (
        genomes_fasta_fork.decompressed.map {
            genome_meta, fasta, index -> [genome_meta, fasta]
        }
    ).mix (
        genomes_pre_processed
    ) .set { ch_genomes }


    SAMTOOLS_FAIDX (
        ch_genomes
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SAM2LCA DB building
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (params.sam2lca_db) {
        sam2lca_db = Channel.fromPath(params.sam2lca_db, checkIfExists: true, type: 'dir')
    } else {
        taxo_nodes = Channel.fromPath(params.taxo_nodes, checkIfExists: true)
        taxo_names = Channel.fromPath(params.taxo_names, checkIfExists: true)
        taxo_merged = Channel.fromPath(params.taxo_merged, checkIfExists: true)
        SAM2LCA_DB (
            ch_genomes,
            taxo_nodes.first(),
            taxo_names.first(),
            taxo_merged.first()
        )

        sam2lca_db = SAM2LCA_DB.out.sam2lca_db
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Mixing read and reference channel
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */


    FASTP.out.reads_merged.mix(FASTP.out.reads) // meta_reads, merged_reads
        .combine(ch_indices) // meta_genome, genome_index
        .map {
            meta_reads, reads, meta_genome, genome_index ->
                [
                    [
                        'id': meta_reads.id + "_" + meta_genome.genome_name,
                        'genome_name': meta_genome.genome_name,
                        'taxid': meta_genome.taxid,
                        'ploidy':meta_genome.ploidy,
                        'sample_name': meta_reads.id,
                        'single_end': params.save_merged ? true : meta_reads.single_end
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

    if (params.estimate_complexity && (! params.deduplicate || params.dedup_tool == 'samtools')) {

        ch_preseq_table = Channel.empty()

        if (params.preseq_mode == 'lc_extrap') {
            PRESEQ_LCEXTRAP(
                ALIGN_BOWTIE2.out.bam
            )

            ch_preseq_table = ch_preseq_table.mix(PRESEQ_LCEXTRAP.out.lc_extrap)

        } else if (params.preseq_mode == 'c_curve') {
            PRESEQ_CCURVE(
                ALIGN_BOWTIE2.out.bam
            )

            ch_preseq_table = ch_preseq_table.mix(PRESEQ_CCURVE.out.c_curve)
        }

        PLOT_PRESEQ(ch_preseq_table)

    }

    if ( params.deduplicate && params.dedup_tool == 'samtools'){
        SAMTOOLS_REMOVE_DUP (
            ALIGN_BOWTIE2.out.bam.join(
                ALIGN_BOWTIE2.out.bai
            )
        )

        bams_synced = SAMTOOLS_REMOVE_DUP.out.bam.map {
            meta, bam -> [['id':meta.sample_name], bam] // meta.id, bam
        }.groupTuple()

    } else {
        bams_synced = ALIGN_BOWTIE2.out.bam.map {
            meta, bam -> [['id':meta.sample_name], bam] // meta.id, bam
        }.groupTuple()
    }

    MERGE_SORT_INDEX_SAMTOOLS (
        bams_synced
    )

    SAM2LCA (
        MERGE_SORT_INDEX_SAMTOOLS.out.bam.join(
            MERGE_SORT_INDEX_SAMTOOLS.out.bai
        ),
        sam2lca_db.first()
    )

    SAM2LCA_MERGE (
        SAM2LCA.out.csv
            .map {
                meta, csv ->
                [ csv ]
            }
            .collect(),
        params.sam2lca_split_rank
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

    BAM_SORT_SAMTOOLS {
        bam_split_by_ref
    }

    BAM_SORT_SAMTOOLS.out.bam.join(
        BAM_SORT_SAMTOOLS.out.bai
    ).map {
        meta, bam, bai -> [meta.taxid, meta, bam, bai] // taxid, meta, bam, bai
    }.combine(
        ch_genomes
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

    synced_ch.dump(tag: 'synced_ch', pretty: true)

    if (params.damage_tool == 'mapdamage2') {
        MAPDAMAGE2 (
            synced_ch.map {
                meta, bam, bai, fasta, fai -> [meta, bam, fasta]
            }
        )
        ch_versions = ch_versions.mix(MAPDAMAGE2.out.versions.first())

        synced_ch = synced_ch.join(
            MAPDAMAGE2.out.rescaled
        ).map {
            meta, bam, bai, fasta, fai, rescaled_bam -> [meta, rescaled_bam, bai, fasta, fai]
        }

    } else {
        DAMAGEPROFILER (
            synced_ch.map {
                meta, bam, bai, fasta, fai -> [meta, bam, fasta, fai]
            },
            []
        )

        ch_versions = ch_versions.mix(DAMAGEPROFILER.out.versions.first())
    }

    QUALIMAP_BAMQC (
        synced_ch
            .map{ meta, bam, bai, fasta, fai  -> [meta, bam] } // meta, bam
    )
    ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions.first())

    if (! params.skip_variant_calling) {
        VARIANT_CALLING (
        synced_ch
    )
        ch_versions = ch_versions.mix(VARIANT_CALLING.out.versions)
    }


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
    if ( ! params.skip_variant_calling) {
        ch_multiqc_files = ch_multiqc_files.mix(VARIANT_CALLING.out.stats.collect{it[1]}.ifEmpty([]))
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
