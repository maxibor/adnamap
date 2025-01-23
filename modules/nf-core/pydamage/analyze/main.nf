process PYDAMAGE_ANALYZE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pysam_numba_pip_pydamage:db4a63d83de6b750' :
        'biocontainers/pydamage:0.90--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("pydamage_results/*_pydamage_results.csv"), emit: csv
    tuple val(meta), path("pydamage_results/plots/*.png"), emit: png, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pydamage \\
        analyze \\
        $args \\
        -p $task.cpus \\
        $bam

    mv pydamage_results/pydamage_results.csv pydamage_results/${prefix}_pydamage_results.csv
    mv pydamage_results/plots/* pydamage_results/plots/${prefix}_pydamage.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pydamage: \$(pydamage --version | sed -n 's/pydamage, version \\(.*\\)/\\1/p')
    END_VERSIONS
    """
}
