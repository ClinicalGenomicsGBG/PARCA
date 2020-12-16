#import glob,re
import yaml
import pandas as pd
import numpy as np
from workflows.utils.process_runinfo_metadata import ProcessRuninfoMetadata

configfile: "config/config.yaml"

singularity: config['singularity_image']

wildcard_constraints:
    start_date="[\d-]+"

def generate_pipeline_input(run_dictionary, out_directory):
    pipeline_input=[]
    for run_id in run_dictionary:
        if run_dictionary[run_id].get("case") and run_dictionary[run_id].get("control"):
            case_sample_id = run_dictionary[run_id].get("case")
            control_sample_id = run_dictionary[run_id].get("control")
            case_control_list=f'{out_directory}/{run_id}/case_control_krona.txt'
            pipeline_input.append(case_control_list)
        elif run_dictionary[run_id].get("case"):
            case_sample_id = run_dictionary[run_id].get("case")
            case_list=f'{out_directory}/{run_id}/case_krona.txt'
            pipeline_input.append(case_list)
    return pipeline_input

run_dict_list = config['run_dict_list']
run_dict = ProcessRuninfoMetadata.nested_run_dict(run_dict_list)

metadata_dict = config['metadata_dict']
metadata_df = pd.DataFrame(metadata_dict)

print(generate_pipeline_input(run_dict, out_directory=config['outdir']))

rule all:
    input:
        #expand("{outdir}/{start_date}_{run_id}/case_control_krona.txt", outdir=config['outdir'], start_date="20201202", run_id="run_1")
        # call the fastqc rule too
        generate_pipeline_input(run_dict, out_directory=config['outdir'])[0]
    run:
        print(input)

rule control_and_case:
    input:
        case = lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage2/kmer_input/kmer_input.fasta".format(
                    sample = ProcessRuninfoMetadata.get_sample(run_dictionary=run_dict,
                                                               run_id=f'{wildcards.start_date}_{wildcards.run_id}',
                                                               case_or_control='case'), 
                    sample_type = ProcessRuninfoMetadata.get_column(df=metadata_df, 
                                                                    run_dictionary=run_dict,
                                                                    run_id=f'{wildcards.start_date}_{wildcards.run_id}',
                                                                    case_or_control="case",
                                                                    column="PE_or_SE",
                                                                    unique=True),
                    nucleotide = ProcessRuninfoMetadata.get_column(df=metadata_df,
                                                                   run_dictionary=run_dict,
                                                                   run_id=f'{wildcards.start_date}_{wildcards.run_id}',
                                                                   case_or_control="case",
                                                                   column="nucleotide",
                                                                   unique=True) ),
        control = lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage2/kmer_input/kmer_input.fasta".format(
                    sample = ProcessRuninfoMetadata.get_sample(run_dictionary=run_dict,
                                                               run_id=f'{wildcards.start_date}_{wildcards.run_id}',
                                                               case_or_control='control'), 
                    sample_type = ProcessRuninfoMetadata.get_column(df=metadata_df,
                                                                    run_dictionary=run_dict,
                                                                    run_id=f'{wildcards.start_date}_{wildcards.run_id}',
                                                                    case_or_control="control",
                                                                    column="PE_or_SE",
                                                                    unique=True),
                    nucleotide = ProcessRuninfoMetadata.get_column(df=metadata_df,
                                                                   run_dictionary=run_dict,
                                                                   run_id=f'{wildcards.start_date}_{wildcards.run_id}',
                                                                   case_or_control="control",
                                                                   column="nucleotide",
                                                                   unique=True) )
    output:
        "{outdir}/{start_date}_{run_id}/case_control_krona.txt"
    shell:
        """
        touch {output}
        """

rule case:
    input:
        case = lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage2/kmer_input/kmer_input.fasta".format(
                    sample = ProcessRuninfoMetadata.get_sample(run_dictionary=run_dict,
                                                               run_id=f'{wildcards.start_date}_{wildcards.run_id}',
                                                               case_or_control='case'), 
                    sample_type = ProcessRuninfoMetadata.get_column(df=metadata_df, 
                                                                    run_dictionary=run_dict,
                                                                    run_id=f'{wildcards.start_date}_{wildcards.run_id}',
                                                                    case_or_control="case",
                                                                    column="PE_or_SE",
                                                                    unique=True),
                    nucleotide = ProcessRuninfoMetadata.get_column(df=metadata_df,
                                                                   run_dictionary=run_dict,
                                                                   run_id=f'{wildcards.start_date}_{wildcards.run_id}',
                                                                   case_or_control="case",
                                                                   column="nucleotide",
                                                                   unique=True) )
    output:
        "{outdir}/{start_date}_{run_id}/case_krona.txt"
    shell:
        """
        touch {output}
        """



# rule all2:
#     input: "{outdir}/201216_test/test.tsv".format(outdir="/apps/bio/dev_repos/parca/demo")


# rule test_create_outdir:
#     output: 
#         "{outdir}/201216_test/test.tsv"
#     log: "{outdir}/201216_test/log/test.log"
#     shell:
#         """
#         mkdir -p {wildcards.outdir}/201216_test
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
