

# nf-core/adnamap pipeline parameters                                                                                                                 
                                                                                                                                                      
aDNA mapping to mulitple genomes                                                                                                                      
                                                                                                                                                      
## Input/output options                                                                                                                               
                                                                                                                                                      
Define where the pipeline should find input data and save output data.                                                                                
                                                                                                                                                      
| Parameter | Description | Type | Default | Required | Hidden |                                                                                      
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                             
| `input` | Path to comma-separated file containing information about the samples in the experiment. <details><summary>Help</summary><small>You will n
| `outdir` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. | `string` | res
| `email` | Email address for completion summary. <details><summary>Help</summary><small>Set this parameter to your e-mail address to get a summary e-
| `multiqc_title` | MultiQC report title. Printed as page header, used for filename if not otherwise specified. | `string` |  |  |  |                 
                                                                                                                                                      
## Reference genomes options                                                                                                                          
                                                                                                                                                      
Reference genomes required for the workflow.                                                                                                          
                                                                                                                                                      
| Parameter | Description | Type | Default | Required | Hidden |                                                                                      
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                             
| `genomes` | Path to genome csv sheet <details><summary>Help</summary><small>A genome sample sheet example can be found in test/data/genomes.csv</sma
                                                                                                                                                      
## Pipeline options                                                                                                                                   
                                                                                                                                                      
Option specific to this pipeline                                                                                                                      
                                                                                                                                                      
| Parameter | Description | Type | Default | Required | Hidden |                                                                                      
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                             
| `save_merged` | Merge forward and reverse file if they overlap for PE data | `string` | True |  |  |                                                
| `save_trimmed_fail` | Also include failed trimmed reads in output | `string` |  |  |  |                                                             
| `deduplicate` | Perform deduplication with fastp (or samtools markdup if complexity estimation is enabled) | `string` | True |  |  |                
| `dedup_tool` | Deduplication tool <details><summary>Help</summary><small>fastp or samtools markdup</small></details>| `string` | fastp |  |  |      
| `estimate_complexity` | Estimate complexity of library using Preseq <details><summary>Help</summary><small>Only when dedup tool in samtools markdup<
| `preseq_mode` | Preseq running mode <details><summary>Help</summary><small>Use c_curve for running on very shallowly sequenced data.<br>Use lc_extra
| `damageprofiler_length` | Window length for DamageProfiler <details><summary>Help</summary><small>Number of bases which are considered for frequency
| `damageprofiler_threshold` | DamageProfiler threshold <details><summary>Help</summary><small>Number of bases which are considered for plotting nucle
                                                                                                                                                      
## sam2lca options                                                                                                                                    
                                                                                                                                                      
Options for Lowest Common Ancestor                                                                                                                    
                                                                                                                                                      
| Parameter | Description | Type | Default | Required | Hidden |                                                                                      
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                             
| `sam2lca_db` | Path to sam2lca database directory | `string` |  |  |  |                                                                             
| `sam2lca_acc2tax` | Type of acc2tax mapping | `string` | nucl |  |  |                                                                               
| `taxo_names` | names.dmp file from NCBI taxonomy | `string` |  |  |  |                                                                              
| `taxo_nodes` | nodes.dmp file from NCBI taxonomy | `string` |  |  |  |                                                                              
| `taxo_merged` | merged.dmp file from NCBI taxonomy | `string` |  |  |  |                                                                            
| `sam2lca_identity` | Sequence identity threshold for LCA | `number` | 0.9 |  |  |                                                                   
| `sam2lca_edit` | Edit distance threshold | `integer` |  |  |  |                                                                                     
| `sam2lca_split_rank` | Rank at which to split bam file per LCA | `string` | species |  |  |                                                         
| `sam2lca_split_n` | Minimum number of reads matching LCA to write alignment file | `integer` | 1000 |  |  |                                         
                                                                                                                                                      
## Variant calling options                                                                                                                            
                                                                                                                                                      
                                                                                                                                                      
                                                                                                                                                      
| Parameter | Description | Type | Default | Required | Hidden |                                                                                      
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                             
| `skip_variant_calling` | Skip variant calling | `boolean` |  |  |  |                                                                                
| `freebayes_min_basequality` | Minimum base quality of base to consider read for variant calling | `integer` | 20 |  |  |                            
| `freebayes_min_mapping_quality` | Minimum mapping quality of base to consider read for variant calling | `integer` | 30 |  |  |                     
| `freebayes_min_coverage` | Minimum coverage at a site to call variant | `integer` | 3 |  |  |                                                       
| `freebayes_minallelefreq` | Minimum alternative allele frequency to call variant | `number` | 0.33 |  |  |                                          
| `freebayes_limit_coverage` | Subsample higher coverage to limit | `integer` | 30 |  |  |                                                            
| `freebayes_report_all_sites` |  | `string` | True |  |  |                                                                                           
| `bcftools_view_high_variant_quality` | High variant quality threshold | `integer` | 30 |  |  |                                                      
| `bcftools_view_medium_variant_quality` | Medium variant quality threshold | `integer` | 20 |  |  |                                                  
| `bcftools_view_minimal_allelesupport` | Minimum allele support to call variant | `integer` | 3 |  |  |                                              
                                                                                                                                                      
## Institutional config options                                                                                                                       
                                                                                                                                                      
Parameters used to describe centralised config profiles. These should not be edited.                                                                  
                                                                                                                                                      
| Parameter | Description | Type | Default | Required | Hidden |                                                                                      
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                             
| `custom_config_version` | Git commit id for Institutional configs. | `string` | master |  | True |                                                  
| `custom_config_base` | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running offline, Nextflow will not
| `config_profile_name` | Institutional config name. | `string` |  |  | True |                                                                        
| `config_profile_description` | Institutional config description. | `string` |  |  | True |                                                          
| `config_profile_contact` | Institutional config contact information. | `string` |  |  | True |                                                      
| `config_profile_url` | Institutional config URL link. | `string` |  |  | True |                                                                     
                                                                                                                                                      
## Max job request options                                                                                                                            
                                                                                                                                                      
Set the top limit for requested resources for any single job.                                                                                         
                                                                                                                                                      
| Parameter | Description | Type | Default | Required | Hidden |                                                                                      
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                             
| `max_cpus` | Maximum number of CPUs that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for t
| `max_memory` | Maximum amount of memory that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit f
| `max_time` | Maximum amount of time that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for t
                                                                                                                                                      
## Generic options                                                                                                                                    
                                                                                                                                                      
Less common options for the pipeline, typically set in a config file.                                                                                 
                                                                                                                                                      
| Parameter | Description | Type | Default | Required | Hidden |                                                                                      
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                             
| `help` | Display help text. | `boolean` |  |  | True |                                                                                              
| `publish_dir_mode` | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The Nextflow `publishDir` opti
| `email_on_fail` | Email address for completion summary, only when pipeline fails. <details><summary>Help</summary><small>An email address to send a 
| `plaintext_email` | Send plain-text email instead of HTML. | `boolean` |  |  | True |                                                               
| `max_multiqc_email_size` | File size limit when attaching MultiQC reports to summary emails. | `string` | 25.MB |  | True |                         
| `monochrome_logs` | Do not use coloured log outputs. | `boolean` |  |  | True |                                                                     
| `multiqc_config` | Custom config file to supply to MultiQC. | `string` |  |  | True |                                                               
| `tracedir` | Directory to keep pipeline Nextflow logs and reports. | `string` | results/pipeline_info |  | True |                                   
| `validate_params` | Boolean whether to validate parameters against the schema at runtime | `boolean` | True |  | True |                             
| `show_hidden_params` | Show all params when using `--help` <details><summary>Help</summary><small>By default, parameters set as _hidden_ in the sche
| `enable_conda` | Run this workflow with Conda. You can also use '-profile conda' instead of providing this parameter. | `boolean` |  |  | True |    
                                                                                                                                                      
## Other parameters                                                                                                                                   
                                                                                                                                                      
| Parameter | Description | Type | Default | Required | Hidden |                                                                                      
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                             
| `schema_ignore_params` |  | `string` | genomes |  |  |                                                                                              
                                                                                                                                                      


