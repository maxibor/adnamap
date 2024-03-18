process SAM2LCA_MERGE {
    label 'process_single'

    conda (params.enable_conda ? "bioconda::sam2lca=1.1.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sam2lca:1.1.4--pyhdfd78af_0' :
        'quay.io/biocontainers/sam2lca:1.1.4--pyhdfd78af_0'            }"

    input:
    path(sam2lca_csv)
    val(rank)

    output:
    path("merged.sam2lca.csv"), emit: sam2lca_summary

    script:
    def args = task.ext.args ?: ""

    """
    merge_sam2lca.py -r $rank
    """
}
