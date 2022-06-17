process SAM2LCA {
    tag "${meta.genome_name}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::sam2lca=1.0.0" : null)
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    //     'ubuntu:20.04' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.sam2lca.json"), emit: json
    tuple val(meta), path("*.sam2lca.csv"),  emit: csv
    tuple val(meta), path("*.sam2lca.bam"),  optional:true, emit: bam
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    gunzip = archive.toString() - '.gz'
    """
    sam2lca \\
        analyze \\
        $args \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sam2lca: \$(echo \$(sam2lca --version 2>&1) | sed 's/sam2lca,\sversion\s//')
    END_VERSIONS
    """
}
