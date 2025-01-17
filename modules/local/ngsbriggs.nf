process NGSBRIGGS {
    tag "${meta.id}"
    label 'process_single'

    conda (params.enable_conda ? "maxibor::ngsbriggs" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/ngsbriggs:1.0.0--42aa832dbea861ba' :
        ''}"

    input:
    tuple val(meta), path(bam), path(fasta)

    output:
    tuple val(meta), path("*_ngsbriggs.txt"), emit: ngsbriggs
    tuple val(meta), path("*_ngsbriggs.log"), emit: ngsbriggs_log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix   = task.ext.prefix ?: "${meta.id}"
    """
    ngsbriggs \\
        -bam $bam \\
        -ref $fasta \\
        $args \\
        -model nb &> ${prefix}_ngsbriggs.log

    grep lambda ${prefix}_ngsbriggs.log > ${prefix}_ngsbriggs.txt
    """
}
