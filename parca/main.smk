#import glob,re
import yaml
import pandas as pd
import numpy as np

configfile: "config/config.yaml"

singularity: config['singularity_image']

def nested_run_dict(list_run_dictionary):
    run_dict_nested = {}
    for run in list_run_dictionary:
        start_date_tmp = run['start_date']
        run_id_tmp = run['run_id']
    
        run_dict_nested[f'{start_date_tmp}_{run_id_tmp}'] = run
    
    return run_dict_nested

def generate_pipeline_input(run_dictionary):
    pipeline_input=[]
    for run_id in run_dictionary:
        if run_dictionary[run_id].get("case") and run_dictionary[run_id].get("control"):
            case_sample_id = run_dictionary[run_id].get("case")
            control_sample_id = run_dictionary[run_id].get("control")
            case_control_list=f'{outdir}/{run_id}/case_control_krona.txt'
            pipeline_input.extend(case_control_list)
        elif run_dictionary[run_id].get("case"):
            case_sample_id = run_dictionary[run_id].get("case")
            case_list=f'{outdir}/{run_id}/case_krona.txt'
            pipeline_input.extend(case_list)
    return pipeline_input

def get_column(df,run_id, case_or_control,column):
    found_column = df.loc[df.get('sample_id') == run_dict.get(run_id).get(case_or_control) ].get(column)
    return found_column

run_dict_list = config['run_dict_list']
run_dict = nested_run_dict(run_dict_list)

metadata_dict = config['metadata_dict']
metadata_df = pd.DataFrame(metadata_dict)

outdir = config['outdir']

rule all:
    input:
        generate_pipeline_input(run_dict)[0]

rule control_and_case:
    input:
        # case = lambda wildcards: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage2/kmer_input/kmer_input.fasta".format(        ,
        #               outdir = wildcards.outdir,
        #               start_date = wildcards.start_date,
        #               run_id = wildcards.run_id,
        #               sample = run_dict[f'{wildcards.start_date}_{wildcards.run_id}']['case'],
        #               sample_type = list(set(metadata_df.loc[(metadata_df['sample_id'] == run_dict[f'{wildcards.start_date}_{wildcards.run_id}']['case'] )]['PE_or_SE']))[0],
        #               nucleotide = list(set(metadata_df.loc[(metadata_df['sample_id'] == run_dict[f'{wildcards.start_date}_{wildcards.run_id}']['case'] )]['nucleotide']))[0]
        #               ),
        # control = expand("{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage2/kmer_input/kmer_input.fasta",
        #               sample = lambda wildcards: run_dict[f'{wildcards.start_date}_{wildcards.run_id}']['control'],
        #               sample_type = lambda wildcards: list(set(metadata_df.loc[(metadata_df['sample_id'] == run_dict[f'{wildcards.start_date}_{wildcards.run_id}']['control'] )]['PE_or_SE']))[0],
        #               nucleotide = lambda wildcards: list(set(metadata_df.loc[(metadata_df['sample_id'] == run_dict[f'{wildcards.start_date}_{wildcards.run_id}']['control'] )]['nucleotide']))[0]
        #               ),
        case = expand("{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage2/kmer_input/kmer_input.fasta",
                      sample = lambda wildcards: run_dict[f'{wildcards.start_date}_{wildcards.run_id}']['case'],
                      sample_type = lambda wildcards: list(set(metadata_df.loc[(metadata_df['sample_id'] == run_dict[f'{wildcards.start_date}_{wildcards.run_id}']['case'] )]['PE_or_SE']))[0],
                      nucleotide = lambda wildcards: list(set(metadata_df.loc[(metadata_df['sample_id'] == run_dict[f'{wildcards.start_date}_{wildcards.run_id}']['case'] )]['nucleotide']))[0]
                      ),
        control = expand("{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage2/kmer_input/kmer_input.fasta",
                      sample = lambda wildcards: run_dict[f'{wildcards.start_date}_{wildcards.run_id}']['control'],
                      sample_type = lambda wildcards: list(set(metadata_df.loc[(metadata_df['sample_id'] == run_dict[f'{wildcards.start_date}_{wildcards.run_id}']['control'] )]['PE_or_SE']))[0],
                      nucleotide = lambda wildcards: list(set(metadata_df.loc[(metadata_df['sample_id'] == run_dict[f'{wildcards.start_date}_{wildcards.run_id}']['control'] )]['nucleotide']))[0]
                      ),
    output:
        "{outdir}/{start_date}_{run_id}/case_control_krona.txt"
    shell:
        """
        touch {output}
        """

# rule case:
#     output:
#         "{outdir}/case_{case}.txt"
#     shell:
#         """
#         touch {output}
#         """

#run_dict = [{'run_id': 'run_1', 'case': 'sample_1', 'control': 'sample_2'}, {'run_id': 'run_2', 'case': 'sample_2'}]
#metadata_dict = [{'sample_id': 'sample_1', 'start_date': 20201104, 'nucleotide': 'RNA', 'fwd_or_rev': 'fwd', 'path_to_file': '/apps/bio/dev_repos/parca/demo/raw_samples/SRR1761912_1.fastq.gz', 'adapters': np.nan, 'PE_or_SE': 'PE'}, {'sample_id': 'sample_1', 'start_date': 20201104, 'nucleotide': 'RNA', 'fwd_or_rev': 'rev', 'path_to_file': '/apps/bio/dev_repos/parca/demo/raw_samples/SRR1761912_2.fastq.gz', 'adapters': np.nan, 'PE_or_SE': 'PE'}, {'sample_id': 'sample_2', 'start_date': 20201104, 'nucleotide': 'DNA', 'fwd_or_rev': 'fwd', 'path_to_file': '/apps/bio/dev_repos/parca/demo/raw_samples/a.fastq.gz', 'adapters': np.nan, 'PE_or_SE': 'SE'}]

# Rule all
# print(expand("{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/all_classed_read_taxid_names.txt",
#             zip,
#             outdir=[config['outdir']]*len(sample_id_list),
#             sample=sample_id_list,
#             sample_type=sample_type_list,
#             nucleotide=nucleotide_list
#             ))


# rule all:
#     input:
#         # expand("{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage2/kmer_input/kmer_input.fasta",
#         #     zip,
#         #     outdir=[config['outdir']]*len(sample_id_list),
#         #     sample= sample_id_list,
#         #     sample_type = sample_type_list,
#         #     nucleotide = nucleotide_list
#         #     )
#         expand("{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/all_classed_read_taxid_names.txt",
#             zip,
#             outdir=[config['outdir']]*len(sample_id_list),
#             sample=sample_id_list,
#             sample_type=sample_type_list,
#             nucleotide=nucleotide_list
#             )

        

##STAGE 1
include:
    "workflows/snakemake_rules/stage1_qc_trim_ec/setup/setup.smk"
include:
    "workflows/snakemake_rules/stage1_qc_trim_ec/quality_control/fastqc.smk"
include:
    "workflows/snakemake_rules/stage1_qc_trim_ec/trimming/bbduk_trimming.smk"
include:
    "workflows/snakemake_rules/stage1_qc_trim_ec/ec_pollux/ec_pollux.smk"
include:
    "workflows/snakemake_rules/stage1_qc_trim_ec/ec_fiona/ec_fiona.smk"

##STAGE 2
include:
    "workflows/snakemake_rules/stage2_assembly/megahit/megahit.smk"
include:
    "workflows/snakemake_rules/stage2_assembly/bbwrap_alignment/bbwrap_alignment.smk"
include:
    "workflows/snakemake_rules/stage2_assembly/merge_contigs_unmapped/merge_contigs_unmapped.smk"

##STAGE 3
include:
    "workflows/snakemake_rules/stage3_kraken_kaiju/kraken_rules/kraken.smk"
include:
    "workflows/snakemake_rules/stage3_kraken_kaiju/kaiju_rules/kaiju.smk"

##STAGE 4
include:
    "workflows/snakemake_rules/stage4_parse_hits/parse_hits.smk"
include:
    "workflows/snakemake_rules/stage4_parse_hits/taxonomy_processing.smk" 

##STAGE 5
include:
    "workflows/snakemake_rules/stage5_blast_processing/blast_processing.smk" 


##STAGE 6
include:
   "workflows/snakemake_rules/stage6_blast_sliced_db/blast_above_species_classed.smk"

##STAGE 7 
include:
   "workflows/snakemake_rules/stage7_blast_remaining_reads/blast_remaining.smk"

##STAGE 8
include:
    "workflows/snakemake_rules/stage8_format_results/format_results.smk"
