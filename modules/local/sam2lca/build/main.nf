process SAM2LCA_BUILD {
    label 'process_single'

    conda (params.enable_conda ? "bioconda::sam2lca=1.1.4--pyhdfd78af_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sam2lca:1.1.4--pyhdfd78af_0' :
        'quay.io/biocontainers/sam2lca:1.1.4--pyhdfd78af_0'            }"

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

    sam2lca -d sam2lca_db \\
        update-db \\
        -t ncbi_local \\
        --taxo_names $taxo_names \\
        --taxo_nodes $taxo_nodes \\
        --taxo_merged $taxo_merged \\
        -a adnamap \\
        --acc2tax_json adnamap.sam2lca.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sam2lca: \$(echo \$(sam2lca --version 2>&1) | sed 's/^sam2lca, version //' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p sam2lca_db
    touch sam2lca_db/test.pkl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sam2lca: \$(echo \$(sam2lca --version 2>&1) | sed 's/^sam2lca, version //' ))
    END_VERSIONS
    """
}
