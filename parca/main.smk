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
            case_control_list=f'{out_directory}/{run_id}/call_case_control_done'
            pipeline_input.append(case_control_list)
        elif run_dictionary[run_id].get("case"):
            case_sample_id = run_dictionary[run_id].get("case")
            case_list=f'{out_directory}/{run_id}/call_case_done'
            pipeline_input.append(case_list)
    return pipeline_input

run_dict_list = config['run_dict_list']
run_dict = ProcessRuninfoMetadata.nested_run_dict(run_dict_list)

metadata_dict = config['metadata_dict']
metadata_df = pd.DataFrame(metadata_dict)

# print(generate_pipeline_input(run_dict, out_directory=config['outdir']))

rule all:
    input:
        #expand("{outdir}/{start_date}_{run_id}/call_case_control_done", outdir=config['outdir'], start_date="20201202", run_id="run_1")
        # call the fastqc rule too
        generate_pipeline_input(run_dict, out_directory=config['outdir'])

rule call_case:
    input: 
        tableview=lambda wildcards: "{webinterface}/{start_date}_{run_id}_web/tableview/case_readcount_tableview.tsv".format(
            webinterface=config['webinterface'],
            start_date=wildcards.start_date,
            run_id=wildcards.run_id
            ),
        krona_html=lambda wildcards: "{webinterface}/{start_date}_{run_id}_web/krona/case.krona.html".format(
            webinterface=config['webinterface'],
            start_date=wildcards.start_date,
            run_id=wildcards.run_id
            ),
        main_page_stats=lambda wildcards: "{webinterface}/{start_date}_{run_id}_web/main_page_stats_case.tsv".format(
            webinterface=config['webinterface'],
            start_date=wildcards.start_date,
            run_id=wildcards.run_id
            ),
        detailed_stats=lambda wildcards: "{webinterface}/{start_date}_{run_id}_web/tableview/case_detailed_stats.tsv".format(
            webinterface=config['webinterface'],
            start_date=wildcards.start_date,
            run_id=wildcards.run_id
            ),
        fastqs=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/fastq_filtering_done".format(
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
        case="{outdir}/{start_date}_{run_id}/call_case_done"
    shell: 
        """
        touch {output.case}
        """

rule call_case_control:
    input: 
        tableview=lambda wildcards: "{webinterface}/{start_date}_{run_id}_web/tableview/case_control_readcount_tableview.tsv".format(
            webinterface=config['webinterface'],
            start_date=wildcards.start_date,
            run_id=wildcards.run_id
            ),
        krona_html=lambda wildcards: "{webinterface}/{start_date}_{run_id}_web/krona/case_control.krona.html".format(
            webinterface=config['webinterface'],
            start_date=wildcards.start_date,
            run_id=wildcards.run_id
            ),
        main_page_stats=lambda wildcards: "{webinterface}/{start_date}_{run_id}_web/main_page_stats_case_control.tsv".format(
            webinterface=config['webinterface'],
            start_date=wildcards.start_date,
            run_id=wildcards.run_id
            ),
        detailed_stats=lambda wildcards: "{webinterface}/{start_date}_{run_id}_web/tableview/case_control_detailed_stats.tsv".format(
            webinterface=config['webinterface'],
            start_date=wildcards.start_date,
            run_id=wildcards.run_id
            ),
        fastqs_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/fastq_filtering_done".format(
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
        fastqs_control=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/fastq_filtering_done".format(
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
        case_control="{outdir}/{start_date}_{run_id}/call_case_control_done"
    shell: 
        """
        touch {output.case_control}
        """

rule link_case_webinterface:
    input: 
        main_page_stats= lambda wildcards:  "{outdir}/{start_date}_{run_id}/main_page_stats_case.tsv".format(
            outdir=config['outdir'],
            start_date=wildcards.start_date,
            run_id=wildcards.run_id 
            ),
        tableview= lambda wildcards:  "{outdir}/{start_date}_{run_id}/tableview/case_readcount_tableview.tsv".format(
            outdir=config['outdir'],
            start_date=wildcards.start_date,
            run_id=wildcards.run_id 
            ),
        krona_html= lambda wildcards: "{outdir}/{start_date}_{run_id}/krona/case.krona.html".format(
            outdir=config['outdir'],
            start_date=wildcards.start_date,
            run_id=wildcards.run_id 
            ),
         detailed_stats_out=lambda wildcards: "{outdir}/{start_date}_{run_id}/tableview/case_detailed_stats.tsv".format(
            outdir=config['outdir'],
            start_date=wildcards.start_date,
            run_id=wildcards.run_id 
            )
    output: 
        main_page_stats= "{webinterface}/{start_date}_{run_id}_web/main_page_stats_case.tsv",
        tableview= "{webinterface}/{start_date}_{run_id}_web/tableview/case_readcount_tableview.tsv",
        krona_html= "{webinterface}/{start_date}_{run_id}_web/krona/case.krona.html",
        detailed_stats="{webinterface}/{start_date}_{run_id}_web/tableview/case_detailed_stats.tsv"
    shell:
        """
        ln -s {input.main_page_stats} {output.main_page_stats}
        ln -s {input.tableview} {output.tableview}
        ln -s {input.krona_html} {output.krona_html}
        ln -s {input.detailed_stats} {output.detailed_stats}
        """

rule link_case_control_webinterface:
    input: 
        main_page_stats= lambda wildcards:  "{outdir}/{start_date}_{run_id}/main_page_stats_case_control.tsv".format(
            outdir=config['outdir'],
            start_date=wildcards.start_date,
            run_id=wildcards.run_id 
            ),
        tableview= lambda wildcards:  "{outdir}/{start_date}_{run_id}/tableview/case_control_readcount_tableview.tsv".format(
            outdir=config['outdir'],
            start_date=wildcards.start_date,
            run_id=wildcards.run_id 
            ),
        krona_html= lambda wildcards: "{outdir}/{start_date}_{run_id}/krona/case_control.krona.html".format(
            outdir=config['outdir'],
            start_date=wildcards.start_date,
            run_id=wildcards.run_id 
            ),
        detailed_stats_out=lambda wildcards: "{outdir}/{start_date}_{run_id}/tableview/case_control_detailed_stats.tsv".format(
            outdir=config['outdir'],
            start_date=wildcards.start_date,
            run_id=wildcards.run_id 
            )
    output: 
        main_page_stats= "{webinterface}/{start_date}_{run_id}_web/main_page_stats_case_control.tsv",
        tableview= "{webinterface}/{start_date}_{run_id}_web/tableview/case_control_readcount_tableview.tsv",
        krona_html= "{webinterface}/{start_date}_{run_id}_web/krona/case_control.krona.html",
        detailed_stats= "{webinterface}/{start_date}_{run_id}_web/tableview/case_control_detailed_stats.tsv"
    shell:
        """
        ln -s {input.main_page_stats} {output.main_page_stats}
        ln -s {input.tableview} {output.tableview}
        ln -s {input.krona_html} {output.krona_html}
        ln -s {input.detailed_stats} {output.detailed_stats}
        """
rule create_main_page_stats_case:
    input: 
        raw_read_count=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage1/samples/count_raw_reads.txt".format(
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
        trimmed_read_count=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage1/trimming/count_bbduk_trimmed_reads.txt".format(
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
        classified_read_count=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage8/tableview/count_classified_reads.txt".format(
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
        unclassified_read_count=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage8/tableview/count_unclassified_reads.txt".format(
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
        stats="{outdir}/{start_date}_{run_id}/main_page_stats_case.tsv"
    shell: 
        """
        echo -e "name\traw_reads\ttrimmed_reads\tclassified_reads\tunclassified_reads" > {output.stats}
        echo -e "{wildcards.start_date}_{wildcards.run_id}\t$(grep "^[0-9]" {input.raw_read_count})\t$(grep "^[0-9]" {input.trimmed_read_count})\t$(grep "^[0-9]" {input.classified_read_count})\t$(grep "^[0-9]" {input.unclassified_read_count})" >> {output.stats}
        """

rule create_main_page_stats_case_control:
    input: 
        case_raw_read_count=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage1/samples/count_raw_reads.txt".format(
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
        case_trimmed_read_count=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage1/trimming/count_bbduk_trimmed_reads.txt".format(
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
        case_classified_read_count=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage8/tableview/count_classified_reads.txt".format(
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
        case_unclassified_read_count=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage8/tableview/count_unclassified_reads.txt".format(
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

        control_raw_read_count=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage1/samples/count_raw_reads.txt".format(
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
                                                                   unique=True) ),
        control_trimmed_read_count=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage1/trimming/count_bbduk_trimmed_reads.txt".format(
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
                                                                   unique=True) ),
        control_classified_read_count=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage8/tableview/count_classified_reads.txt".format(
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
                                                                   unique=True) ),
        control_unclassified_read_count=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage8/tableview/count_unclassified_reads.txt".format(
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
        stats="{outdir}/{start_date}_{run_id}/main_page_stats_case_control.tsv"
    shell: 
        """
        raw_read_count=$(echo $(grep "^[0-9]" {input.case_raw_read_count})+$(grep "^[0-9]" {input.control_raw_read_count})|bc)
        trimmed_read_count=$(echo $(grep "^[0-9]" {input.case_trimmed_read_count})+$(grep "^[0-9]" {input.control_trimmed_read_count})|bc)
        classified_read_count=$(echo $(grep "^[0-9]" {input.case_classified_read_count})+$(grep "^[0-9]" {input.control_classified_read_count})|bc)
        unclassified_read_count=$(echo $(grep "^[0-9]" {input.case_unclassified_read_count})+$(grep "^[0-9]" {input.control_unclassified_read_count})|bc)
        
        echo -e "name\traw_reads\ttrimmed_reads\tclassified_reads\tunclassified_reads" > {output.stats}
        echo -e "{wildcards.start_date}_{wildcards.run_id}\t$raw_read_count\t$trimmed_read_count\t$classified_read_count\t$unclassified_read_count" >> {output.stats}
        """

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