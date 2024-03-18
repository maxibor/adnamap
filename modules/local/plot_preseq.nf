process PLOT_PRESEQ {
    tag "${meta.id}-${meta.genome_name}"

    conda (params.enable_conda ? "conda-forge:matplotlib=3.5.2 conda-forge:pandas=1.4.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d0aaa59a9c102cbb5fb1b38822949827c5119e45:a53fc5f5fda400a196436eac5c44ff3e2d42b0dc-0' :
        'quay.io/biocontainers/mulled-v2-d0aaa59a9c102cbb5fb1b38822949827c5119e45:a53fc5f5fda400a196436eac5c44ff3e2d42b0dc-0'            }"

    input:
    tuple val(meta), path(preseq_table)

    output:
    path '*.png'       , emit: png

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def args2 = task.ext.args2 ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    script: // This script is bundled with the pipeline, in nf-core/adnamap/bin/
    """
    plot_preseq.py \\
        $args \\
        -s $prefix \\
        -o ${prefix}.png \\
        $preseq_table
    """
}
