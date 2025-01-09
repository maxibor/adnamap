//
// Check input samplesheet and get read channels
//

include { GUNZIP        } from '../../modules/nf-core/gunzip/main'
include { UNTAR         } from '../../modules/nf-core/untar/main'
include { BOWTIE2_BUILD } from '../../modules/nf-core/bowtie2/build/main'

workflow GENOME_CHECK {
    take:
    genome_sheet 

    main:

    genome_sheet
        .splitCsv(header:true, sep:',')
        .map { row -> create_genome_channel(row) }
        .set { genomes }

    genomes
        .branch {
            decompressed: it[1].getExtension() != 'gz' 
            compressed: it[1].getExtension() == 'gz'
        }
        .set { genomes_fasta_fork }

    GUNZIP (
        genomes_fasta_fork.compressed
        .map {
            genome_meta, fasta, index -> [genome_meta, fasta]
        }
    ) // decompress the genome fasta

    GUNZIP.out.gunzip
        .mix(
                genomes_fasta_fork.decompressed.map {
                    genome_meta, fasta, index -> [genome_meta, fasta]
                }
        ).set { genomes_pre_processed } // genome_fasta decompressed

    genomes_pre_processed.join(
        genomes
    ).map {
        genome_meta, fasta, fasta_raw, index_raw -> [genome_meta, fasta, index_raw] 
    }.branch {
            no_index: ! it[2] // decompressed genome, no idx
            has_index_decompressed: it[2] && it[2].getExtension() != 'gz' // genome with decompressed index
            has_index_compressed: it[2] && it[2].getExtension() == 'gz' // genome with compressed index
    } .set { genomes_idx_fork }

    UNTAR(
        genomes_idx_fork.has_index_compressed.map {
            meta, fasta, index_raw ->
            [meta, index_raw]
        }
    )

    BOWTIE2_BUILD (
        genomes_idx_fork.no_index.map {
            meta, fasta, index_raw ->
            [meta, fasta]
        }
    )

    ch_indices = BOWTIE2_BUILD.out.index.mix(
            genomes_idx_fork.has_index_decompressed
            .map {
                meta_genome, genome_fasta, genome_index ->
                [meta_genome, genome_index]
            }
        ).mix (
                UNTAR.out.untar
        ) // all genomes indices ready
    

    emit:
    genomes_pre_processed                                // channel: [ [genome_meta], [ genome_fasta ] ]
    ch_indices                                // channel: [ [genome_meta], [ genome_index ] ]
}

def create_genome_channel(LinkedHashMap row) {
    def meta = [:]
    meta.genome_name = row.genome_name
    meta.taxid       = row.taxid
    meta.ploidy      = row.ploidy

    // add path(s) of the fastq file(s) to the meta map
    def genome_meta = []
    if (!file(row.genome_path).exists()) {
        exit 1, "ERROR: Please check input genome sheet -> Genome path does not exist!\n${row.genome_path}"
    }
    if (row.genome_index != "") {
       if (file(row.genome_index).exists() and ! file(row.genome_index).isEmpty()) {
            genome_meta = [ meta, file(row.genome_path), file(row.genome_index) ]
       } else {
            genome_meta = [ meta, file(row.genome_path), null ]
       }
    } else {
        genome_meta = [ meta, file(row.genome_path), null ]
    }
    return genome_meta
}
