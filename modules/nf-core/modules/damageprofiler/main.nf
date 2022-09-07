process DAMAGEPROFILER {
    tag "${meta.id}"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::damageprofiler=1.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/damageprofiler:1.1--hdfd78af_2' :
        'quay.io/biocontainers/damageprofiler:1.1--hdfd78af_2' }"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta), path(fai)

    output:
    tuple val(meta), path("${prefix}_damageprofiler"), emit: results
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mv $bam ${prefix}.sam2lca.bam
    mv $bai ${prefix}.sam2lca.bai

    damageprofiler \\
        -Xmx${task.memory.toMega()}m \\
        -i ${prefix}.sam2lca.bam \\
        -r $fasta \\
        -title $prefix \\
        -l ${params.damageprofiler_length} \\
        -t ${params.damageprofiler_threshold} \\
        -o ${prefix}_damageprofiler \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        damageprofiler: \$(damageprofiler -v | sed 's/^DamageProfiler v//')
    END_VERSIONS
    """
}
