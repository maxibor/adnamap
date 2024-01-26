include { CREATE_ACC2TAX } from '../../modules/local/create_acc2tax'
include { SAM2LCA_BUILD } from '../../modules/local/sam2lca/build/main'

workflow SAM2LCA_DB {
    take:
        genomes // meta, fasta
        taxo_nodes // nodes.dmp
        taxo_names // names.dmp
        taxo_merged // merged.dmp

    main:
        CREATE_ACC2TAX(genomes)

        acc2tax = CREATE_ACC2TAX.out.acc2tax.collectFile(
            name: 'adnamap.accession2taxid',
            keepHeader: true
        )

        SAM2LCA_BUILD (
            acc2tax,
            taxo_nodes,
            taxo_names,
            taxo_merged
        )

    emit:
        sam2lca_db = SAM2LCA_BUILD.out.sam2lca_db
}
