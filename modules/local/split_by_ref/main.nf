process SPLIT_BY_REF {

    tag "${meta.genome_name}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::pysam=0.19.1" : null)
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    //     'ubuntu:20.04' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(met), path("*.bam"), emit: bam

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python \\
        split_bam_by_taxid.py \\
        -t ${task.cpus} \\
        $bam
    """
}
