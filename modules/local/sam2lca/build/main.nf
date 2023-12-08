process SAM2LCA_BUILD {
    label 'process_single'

    conda (params.enable_conda ? "bioconda::sam2lca=1.1.3=pyhdfd78af_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sam2lca:1.1.3--pyhdfd78af_0' :
        'quay.io/biocontainers/sam2lca:1.1.3--pyhdfd78af_0'            }"

    input:
    path(acc2tax)
    path(taxo_nodes) // nodes.dmp
    path(taxo_names) // names.dmp
    path(taxo_merged) // merged.dmp

    output:
    path("sam2lca_db"), emit: sam2lca_db

    script:
    def args = task.ext.args ?: ""

    """
    mkdir -p sam2lca_db
    gzip $acc2tax
    md5sum ${acc2tax}.gz > ${acc2tax}.gz.md5
    sam2lca_json.py ${acc2tax}.gz ${acc2tax}.gz.md5

    sam2lca -t ncbi_local \\
        --taxo_names $taxo_names \\
        --taxo_nodes $taxo_nodes \\
        --taxo_merged $taxo_merged \\
        -a adnamap \\
        --acc2tax_json adnamap.sam2lca.json
    """
}
