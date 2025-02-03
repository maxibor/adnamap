process COMPUTE_LAMBDA {
    tag "${meta.id}"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::sam2lca=1.1.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_matplotlib_numpy_pandas_scipy:1777c3fa3e3667ed' :
        'quay.io/biocontainers/sam2lca:1.1.4--pyhdfd78af_0'            }"

    input:
    tuple val(meta), path(length)

    output:
    tuple val(meta), path("*.tsv"), emit: lambda_tsv
    tuple val(meta), path("*.png"), emit: png

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    compute_lambda_kistler.py \\
        $length \\
        $prefix
    """
}
