process BOWTIE2_BUILD {
    tag "$meta.genome_name"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::bowtie2=2.4.4' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.4.4--py39hbb4e92a_0' :
        'quay.io/biocontainers/bowtie2:2.4.4--py39hbb4e92a_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('bowtie2'), emit: index
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def basename = fasta.baseName.toString().tokenize(".")[0]
    """
    mkdir bowtie2
    bowtie2-build $args --threads $task.cpus $fasta bowtie2/$basename
    mv $fasta bowtie2
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS
    """
}
