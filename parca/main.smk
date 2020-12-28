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
        generate_pipeline_input(run_dict, out_directory=config['outdir'])

# rule control_and_case:
#     input:
#         case = lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/readcount_tableview.tsv".format(
#                     sample = ProcessRuninfoMetadata.get_sample(run_dictionary=run_dict,
#                                                                run_id=f'{wildcards.start_date}_{wildcards.run_id}',
#                                                                case_or_control='case'), 
#                     sample_type = ProcessRuninfoMetadata.get_column(df=metadata_df, 
#                                                                     run_dictionary=run_dict,
#                                                                     run_id=f'{wildcards.start_date}_{wildcards.run_id}',
#                                                                     case_or_control="case",
#                                                                     column="PE_or_SE",
#                                                                     unique=True),
#                     nucleotide = ProcessRuninfoMetadata.get_column(df=metadata_df,
#                                                                    run_dictionary=run_dict,
#                                                                    run_id=f'{wildcards.start_date}_{wildcards.run_id}',
#                                                                    case_or_control="case",
#                                                                    column="nucleotide",
#                                                                    unique=True) ),
#         control = lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/organism_dfs".format(
#                     sample = ProcessRuninfoMetadata.get_sample(run_dictionary=run_dict,
#                                                                run_id=f'{wildcards.start_date}_{wildcards.run_id}',
#                                                                case_or_control='control'), 
#                     sample_type = ProcessRuninfoMetadata.get_column(df=metadata_df,
#                                                                     run_dictionary=run_dict,
#                                                                     run_id=f'{wildcards.start_date}_{wildcards.run_id}',
#                                                                     case_or_control="control",
#                                                                     column="PE_or_SE",
#                                                                     unique=True),
#                     nucleotide = ProcessRuninfoMetadata.get_column(df=metadata_df,
#                                                                    run_dictionary=run_dict,
#                                                                    run_id=f'{wildcards.start_date}_{wildcards.run_id}',
#                                                                    case_or_control="control",
#                                                                    column="nucleotide",
#                                                                    unique=True) )
#     output:
#         "{outdir}/{start_date}_{run_id}/case_control_krona.txt"
#     shell:
#         """
#         touch {output}
#         """

# rule case:
#     input:
#         case = lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/all_classed_read_taxid_names.txt".format(
#                     sample = ProcessRuninfoMetadata.get_sample(run_dictionary=run_dict,
#                                                                run_id=f'{wildcards.start_date}_{wildcards.run_id}',
#                                                                case_or_control='case'), 
#                     sample_type = ProcessRuninfoMetadata.get_column(df=metadata_df, 
#                                                                     run_dictionary=run_dict,
#                                                                     run_id=f'{wildcards.start_date}_{wildcards.run_id}',
#                                                                     case_or_control="case",
#                                                                     column="PE_or_SE",
#                                                                     unique=True),
#                     nucleotide = ProcessRuninfoMetadata.get_column(df=metadata_df,
#                                                                    run_dictionary=run_dict,
#                                                                    run_id=f'{wildcards.start_date}_{wildcards.run_id}',
#                                                                    case_or_control="case",
#                                                                    column="nucleotide",
#                                                                    unique=True) )
#     output:
#         "{outdir}/{start_date}_{run_id}/case_krona.txt"
#     shell:
#         """
#         touch {output}
#         """

rule call_case:
    input: 
        tableview="{outdir}/{start_date}_{run_id}/tableview/case_readcount_tableview.tsv",
        krona_html="{outdir}/{start_date}_{run_id}/krona/case.krona.html"
    output: ""
    shell: ""


rule call_case_control:
    input: 
        tableview="{outdir}/{start_date}_{run_id}/tableview/case_control_readcount_tableview.tsv",
        krona_html="{outdir}/{start_date}_{run_id}/krona/case_control.krona.html"
    output: ""
    shell: ""



rule generate_krona_plot_case:
    input:
        readcount_krona = lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/krona/readcount_krona.tsv".format(
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
        krona_html="{outdir}/{start_date}_{run_id}/krona/case.krona.html"
    conda: "../../conda/krona.yaml"
    shell:
        """
        ktImportText {input.readcount_krona},"case" -o {output.krona_html};
        """

rule generate_krona_plot_case_control:
    input:
        readcount_krona_case = lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/krona/readcount_krona.tsv".format(
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
        readcount_krona_control = lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/krona/readcount_krona.tsv".format(
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
        krona_html="{outdir}/{start_date}_{run_id}/krona/case_control.krona.html"
    conda: "../../conda/krona.yaml"
    shell:
        """
        ktImportText {input.readcount_krona_case},"case" {input.readcount_krona_control},"control" -o {output.krona_html};
        """

rule tableview_case_control:
    input: 
        case = lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/{sample_type}_{nucleotide}stage8/tableview/readcount_tableview.tsv".format(
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
        control = lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/readcount_tableview.tsv".format(
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
        case_control="{outdir}/{start_date}_{run_id}/tableview/case_control_readcount_tableview.tsv"
    conda: "../../conda/R_env.yaml"
    script: "../../scripts/reformat_results/tableview_case_control.R"   

rule tableview_case:
    input: 
        case="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/readcount_tableview.tsv"
    output: 
        case="{outdir}/{start_date}_{run_id}/tableview/case_readcount_tableview.tsv"
    shell: 
        """
        cp {input.case} {output.case}
        """  

rule name:
    input: 
        raw_read_count="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage1/samples/count_raw_reads.txt",
        trimmed_read_count="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage1/trimming/count_bbduk_trimmed_reads.txt",
        classified_read_count="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage8/tableview/count_classified_reads.txt",
        unclassified_read_count="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage8/tableview/count_unclassified_reads.txt",
        fastq_fitering="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}stage8/tableview/fastq_filtering_done"
    output: "{outdir}/{start_date}_{run_id}"
    run: 

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
include:
    "workflows/snakemake_rules/stage8_format_results/krona_plot.smk"
include:
    "workflows/snakemake_rules/stage8_format_results/generate_files_for_download.smk"