process SAM2LCA {
    tag "${meta.id}"

    conda (params.enable_conda ? "bioconda::sam2lca=1.1.4--pyhdfd78af_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sam2lca:1.1.4--pyhdfd78af_0' :
        'quay.io/biocontainers/sam2lca:1.1.4--pyhdfd78af_0'            }"

    input:
    tuple val(meta), path(bam), path(bai)
    path sam2lca_db

    output:
    tuple val(meta), path("*.sam2lca.json"), emit: json
    tuple val(meta), path("*.sam2lca.csv"),  emit: csv
    tuple val(meta), path("*.sam2lca.bam"),  optional:true, emit: bam
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def args2 = task.ext.args2 ?: ""
    def args3 = task.ext.args3 ?: ""
    """
    sam2lca \\
        -d $sam2lca_db \\
        analyze \\
        -p ${task.cpus} \\
        -b \\
        $args \\
        $args2 \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sam2lca: \$(echo \$(sam2lca --version 2>&1) | sed 's/sam2lca,\sversion\s//')
    END_VERSIONS
    """
}
