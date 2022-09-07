# nf-core/adnamap pipeline parameters

aDNA mapping to mulitple genomes

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `input` | Path to comma-separated file containing information about the samples in the experiment. 
<details><summary>Help</summary><small>You will need to create a design file with information about the samples in your experiment 
before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a 
header row. See (https://nf-co.re/adnamap/usage#samplesheet-input).</small></details>| `string` |  |  |  |
| `outdir` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud 
infrastructure. | `string` |  |  |  |
| `email` | Email address for completion summary. <details><summary>Help</summary><small>Set this parameter to your e-mail address 
to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file 
(`~/.nextflow/config`) then you don't need to specify this on the command line for every run.</small></details>| `string` |  |  |  
|
| `multiqc_title` | MultiQC report title. Printed as page header, used for filename if not otherwise specified. | `string` |  |  | 
|

## Reference genomes options

Reference genomes required for the workflow.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `genomes` | Path to genome csv sheet <details><summary>Help</summary><small>A genome sample sheet example can be found in 
test/data/genomes.csv</small></details>| `string` |  |  |  |

## Pipeline options

Option specific to this pipeline

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `save_merged` | Merge forward and reverse file if they overlap for PE data | `string` | true |  |  |
| `save_trimmed_fail` | Also include failed trimmed reads in output | `string` | false |  |  |
| `deduplicate` | Perform deduplication with fastp | `string` | true |  |  |
| `damageprofiler_length` | Window length for DamageProfiler <details><summary>Help</summary><small>Number of bases which are 
considered for frequency computations</small></details>| `integer` | 100 |  |  |
| `damageprofiler_threshold` | DamageProfiler threshold <details><summary>Help</summary><small>Number of bases which are considered
for plotting nucleotide misincorporations</small></details>| `integer` | 15 |  |  |

## sam2lca options

Options for Lowest Common Ancestor

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `sam2lca_db` | Path to sam2lca database directory | `string` | /Users/maxime/.sam2lca |  |  |
| `sam2lca_acc2tax` | Type of acc2tax mapping | `string` | nucl |  |  |
| `sam2lca_identity` |  | `number` | 0.9 |  |  |
| `sam2lca_split_rank` | Rank at which to split bam file per LCA | `string` | species |  |  |
| `sam2lca_split_n` | Minimum number of reads matching LCA to write alignment file | `integer` | 1000 |  |  |

## Variant calling options



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `freebayes_ploidy` | Ploidy of the reference organism | `integer` | 1 |  |  |
| `freebayes_min_basequality` | Minimum base quality of base to consider read for variant calling | `integer` | 20 |  |  |
| `freebayes_min_mapping_quality` | Minimum mapping quality of base to consider read for variant calling | `integer` | 30 |  |  |
| `freebayes_min_coverage` | Minimum coverage at a site to call variant | `integer` | 3 |  |  |
| `freebayes_minallelefreq` | Minimum alternative allele frequency to call variant | `number` | 0.33 |  |  |
| `freebayes_limit_coverage` | Subsample higher coverage to limit | `integer` | 30 |  |  |
| `freebayes_report_all_sites` |  | `string` | true |  |  |
| `bcftools_view_high_variant_quality` | High variant quality threshold | `integer` | 30 |  |  |
| `bcftools_view_medium_variant_quality` | Medium variant quality threshold | `integer` | 20 |  |  |
| `bcftools_view_minimal_allelesupport` | Minimum allele support to call variant | `integer` | 3 |  |  |

## Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `custom_config_version` | Git commit id for Institutional configs. | `string` | master |  | True |
| `custom_config_base` | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running 
offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't need them, then this is 
not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this 
parameter.</small></details>| `string` | https://raw.githubusercontent.com/nf-core/configs/master |  | True |
| `config_profile_name` | Institutional config name. | `string` |  |  | True |
| `config_profile_description` | Institutional config description. | `string` |  |  | True |
| `config_profile_contact` | Institutional config contact information. | `string` |  |  | True |
| `config_profile_url` | Institutional config URL link. | `string` |  |  | True |

## Max job request options

Set the top limit for requested resources for any single job.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `max_cpus` | Maximum number of CPUs that can be requested for any single job. <details><summary>Help</summary><small>Use to set 
an upper-limit for the CPU requirement for each process. Should be an integer e.g. `--max_cpus 1`</small></details>| `integer` | 16
|  | True |
| `max_memory` | Maximum amount of memory that can be requested for any single job. <details><summary>Help</summary><small>Use to 
set an upper-limit for the memory requirement for each process. Should be a string in the format integer-unit e.g. `--max_memory 
'8.GB'`</small></details>| `string` | 128.GB |  | True |
| `max_time` | Maximum amount of time that can be requested for any single job. <details><summary>Help</summary><small>Use to set 
an upper-limit for the time requirement for each process. Should be a string in the format integer-unit e.g. `--max_time 
'2.h'`</small></details>| `string` | 240.h |  | True |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `help` | Display help text. | `boolean` |  |  | True |
| `publish_dir_mode` | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The 
Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the 
pipeline what method should be used to move these files. See [Nextflow 
docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>| `string` | copy |  | True |
| `email_on_fail` | Email address for completion summary, only when pipeline fails. <details><summary>Help</summary><small>An email
address to send a summary email to when the pipeline is completed - ONLY sent if the pipeline does not exit 
successfully.</small></details>| `string` |  |  | True |
| `plaintext_email` | Send plain-text email instead of HTML. | `boolean` |  |  | True |
| `max_multiqc_email_size` | File size limit when attaching MultiQC reports to summary emails. | `string` | 25.MB |  | True |
| `monochrome_logs` | Do not use coloured log outputs. | `boolean` |  |  | True |
| `multiqc_config` | Custom config file to supply to MultiQC. | `string` |  |  | True |
| `tracedir` | Directory to keep pipeline Nextflow logs and reports. | `string` | ${params.outdir}/pipeline_info |  | True |
| `validate_params` | Boolean whether to validate parameters against the schema at runtime | `boolean` | True |  | True |
| `show_hidden_params` | Show all params when using `--help` <details><summary>Help</summary><small>By default, parameters set as 
_hidden_ in the schema are not shown on the command line when a user runs with `--help`. Specifying this option will tell the 
pipeline to show all parameters.</small></details>| `boolean` |  |  | True |
| `enable_conda` | Run this workflow with Conda. You can also use '-profile conda' instead of providing this parameter. | `boolean`
|  |  | True |

