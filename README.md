# aDNAMap

**adnamap** is a bioinformatics best-practice analysis pipeline for aDNA mapping to mulitple genomes.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/adnamap/results).

## Pipeline summary

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```console
   nextflow run nf-core/adnamap -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!


   ```console
   nextflow run nf-core/adnamap --input samplesheet.csv  --genomes genome_sheet.csv --outdir <OUTDIR> -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
   ```

## Documentation

The adnamap pipeline comes with documentation about the pipeline [usage](docs/usage), [parameters](docs/parameters) and [output](docs/output).

## Credits

nf-core/adnamap was originally written by Maxime Borry.
## Pipeline workflow

```mermaid
flowchart TD
    subgraph fastq_process[FastQ preprocessing]
        fastq{FastQ file}
        fastqc_before[FastQC before]
        fastp["fastp: \n Adapter and Quality trimming + merging + deduplication"]
        fastqc_after[FastQC after]
        fastq --> fastqc_before
        fastq --> fastp
        fastp -- Trimmed fastq --> fastqc_after
    end
    subgraph fastq_preproprecessing[FastA preprocessing]
        fasta{Genome FastA file}
        faidx[Samtools faidx]
        fasta --compressed-->gunzip
        gunzip-->bowtie2_build
        bowtie2_build[Bowtie2-build]
        fasta--not compressed--> bowtie2_build
        gunzip-->faidx
        fasta--not compressed-->faidx
    end
    subgraph alignment[Alignment]
        bowtie2_align["Bowtie2 align subworkflow"]
        bowtie2_build --Bowtie 2 index--> bowtie2_align
        fastp --Trimmed fastq--> bowtie2_align
    end
    subgraph lca[LCA]
      bowtie2_align --samtools merge--> merged_BAM["Merged BAM per sample"]
      merged_BAM --sam2lca--> rank_BAM["taxon rank specific BAM"]
    end
    subgraph alignment_post[Alignment post_processing]
        qualimap[Alignment stats reporting]
        damageprofiler[DamageProfiler: aDNA stats]
        fasta--decompressed-->damageprofiler
        gunzip-->damageprofiler
        rank_BAM--BAM+index-->damageprofiler
        rank_BAM--BAM+index-->qualimap
    end
    subgraph variant_calling[Variant calling]
        snpAD["snpAD: ancient DNA damage aware genotyper"]
        freebayes["Freebayes: genotyper"]

        rank_BAM--BAM+index-->snpAD
        rank_BAM--BAM+index-->freebayes
        fasta--decompressed-->snpAD
        gunzip-->freebayes
        faidx-->freebayes
        fasta--decompressed-->snpAD
        gunzip-->freebayes
        faidx-->snpAD
        vcf{VCF}
        freebayes --> vcf
        vcf-- bcftools consensus -->cons_fa["Consensus genome fasta"]
    end
```

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/adnamap for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline integrates code templates and was made with tools from the [`nf-core` community](https://nf-co.re/). More in their publication below:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
