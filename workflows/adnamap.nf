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
include { BAM_SORT_SAMTOOLS as BSS  } from '../subworkflows/nf-core/bam_sort_samtools/main'


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
include { QUALIMAP_BAMQC                                   } from '../modules/nf-core/modules/qualimap/bamqc/main'
include { DAMAGEPROFILER                                   } from '../modules/nf-core/modules/damageprofiler/main'
include { FREEBAYES                                        } from '../modules/nf-core/modules/freebayes/main'
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
                [meta_reads + meta_genome, reads, genome_index]
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
        it -> [['id':it[0].id], it[1]] // id, bam
    }.groupTuple().set { bams_synced }

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
            [new_meta , bam]
    }
    .set{ bam_split_by_ref} // id, taxid, bam


    INDEX_PER_GENOME {
        bam_split_by_ref
    }

    bam_split_by_ref.join(
        INDEX_PER_GENOME.out.bai
    ).map {
        it -> [it[0].taxid, it[0].id, it[1], it[2]] // taxid, id, bam, bai
    }.combine(
        genomes_pre_processed
            .map{
                it -> [it[0].taxid, it[1]] //taxid, fasta
            }
    , by: 0).combine(                      // taxid, id, bam, bai, fasta
        SAMTOOLS_FAIDX.out.fai
            .map{
                it -> [it[0].taxid, it[0].genome_name, it[1]] // taxid, genome_name, fai
            }
    , by: 0).map{ //taxid, id, bam, bai, fasta, genome_name, fai
        it -> [['id':it[1], 'taxid':it[0], 'genome_name':it[5]], it[2], it[3], it[4], it[6]] // meta, bam, bai, fasta, fai
    }.set {
        synced_ch
    }


    DAMAGEPROFILER (
        synced_ch
    )

    ch_versions = ch_versions.mix(DAMAGEPROFILER.out.versions.first())

    FREEBAYES (
        synced_ch
    )

    ch_versions = ch_versions.mix(FREEBAYES.out.versions.first())

    QUALIMAP_BAMQC (
        synced_ch
            .map{ it -> [it[0], it[1]] } // meta, bam
    )
    ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions.first())




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
