# from ...utils.process_runinfo_metadata import ProcessRuninfoMetadata
import sys
sys.path.append("../..") 

checkpoint tableview:
    input: 
        read_count="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/readcount.tsv",
        trimmed_read_count="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage1/trimming/count_bbduk_trimmed_reads.txt"
    output: 
        tableview="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/readcount_tableview.tsv",
        classified_reads_mincount="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/classified_reads_mincount.tsv",
        read_count="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage8/tableview/count_classified_reads.txt",
        organism_dir=directory("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/organism_dir"),
        kingdom_dir=directory("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/kingdom_dir")
    params:
        # SE_or_PE="SE",
        mincount=config['tableview_min_count']
    conda: "../../conda/R_env.yaml"
    script: "../../scripts/reformat_results/tableview_splitting.R"

rule filter_fastq_organism_SE:
    input: 
        organism_tableview="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage8/tableview/organism_dir/organism_{taxid}.tsv",
        trimmed_reads="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage1/trimming/trimmed_reads.fq"
    output: 
        fastq_out="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage8/tableview/organism_fastq/{taxid}.fastq"
    params:
        SE_or_PE="SE",
        negate_query=False
    conda: "../../conda/R_env.yaml"
    script: "../../scripts/reformat_results/filter_fastq.R"

rule filter_fastq_organism_PE:
    input: 
        organism_tableview="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage8/tableview/organism_dir/organism_{taxid}.tsv",
        trimmed_reads_unmerged="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/unmerged_reads_trimmed.fq",
        trimmed_reads_merged="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/merged_reads_trimmed.fq"
    output: 
        fastq_out="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage8/tableview/organism_fastq/{taxid}.fastq"
    params:
        SE_or_PE="PE",
        negate_query=False
    conda: "../../conda/R_env.yaml"
    script: "../../scripts/reformat_results/filter_fastq.R"


rule filter_fastq_kingdom_SE:
    input: 
        organism_tableview="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage8/tableview/kingdom_dir/kingdom_{kingdom}.tsv",
        trimmed_reads="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage1/trimming/trimmed_reads.fq"
    output: 
        fastq_out="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage8/tableview/kingdom_fastq/{kingdom}.fastq"
    params:
        SE_or_PE="SE",
        negate_query=False
    conda: "../../conda/R_env.yaml"
    script: "../../scripts/reformat_results/filter_fastq.R"

rule filter_fastq_kingdom_PE:
    input: 
        organism_tableview="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage8/tableview/kingdom_dir/kingdom_{kingdom}.tsv",
        trimmed_reads_unmerged="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/unmerged_reads_trimmed.fq",
        trimmed_reads_merged="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/merged_reads_trimmed.fq"
    output: 
        fastq_out="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage8/tableview/kingdom_fastq/{kingdom}.fastq"
    params:
        SE_or_PE="PE",
        negate_query=False
    conda: "../../conda/R_env.yaml"
    script: "../../scripts/reformat_results/filter_fastq.R"

rule filter_fastq_unclassified_SE:
    input: 
        organism_tableview="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage8/tableview/classified_reads_mincount.tsv",
        trimmed_reads="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage1/trimming/trimmed_reads.fq"
    output: 
        fastq_out="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage8/tableview/unclassified_fastq/unclassified.fastq"
    params:
        SE_or_PE="SE",
        negate_query=True
    conda: "../../conda/R_env.yaml"
    threads: 10
    script: "../../scripts/reformat_results/filter_fastq.R"

rule filter_fastq_unclassified_PE:
    input: 
        organism_tableview="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage8/tableview/classified_reads_mincount.tsv",
        trimmed_reads_unmerged="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/unmerged_reads_trimmed.fq",
        trimmed_reads_merged="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/merged_reads_trimmed.fq"
    output: 
        fastq_out="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage8/tableview/unclassified_fastq/unclassified.fastq"
    params:
        SE_or_PE="PE",
        negate_query=True
    conda: "../../conda/R_env.yaml"
    threads: 8
    script: "../../scripts/reformat_results/filter_fastq.R"

rule zip_filtered_fastq_organism:
    input: 
        fastq="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/organism_fastq/{taxid}.fastq"
    output: 
        fastq="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/organism_fastq/{taxid}.fastq.gz"
    threads: 4
    shell:
        """
        pigz -p {threads} -k {input.fastq};
        """ 

rule zip_filtered_fastq_kingdom:
    input: 
        fastq="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/kingdom_fastq/{kingdom}.fastq"
    output: 
        fastq="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/kingdom_fastq/{kingdom}.fastq.gz"
    threads: 4
    shell:
        """
        pigz -p {threads} -k {input.fastq};
        """ 

rule zip_filtered_fastq_unclassified:
    input: 
        fastq="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/unclassified_fastq/unclassified.fastq"
    output: 
        fastq="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/unclassified_fastq/unclassified.fastq.gz",
        read_count="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage8/tableview/count_unclassified_reads.txt"
    threads: 8
    shell:
        """
        echo count > {output.read_count};
        echo $(cat {input.fastq}|wc -l)/4|bc  >> {output.read_count};
        pigz -p {threads} -k {input.fastq};
        """ 

rule link_filtered_fastq_organism:
    input: 
        fastq= lambda wildcards: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/organism_fastq/{taxid}.fastq.gz".format(
            outdir=config['outdir'],
            start_date=wildcards.start_date,
            run_id=wildcards.run_id,
            sample=wildcards.sample,
            sample_type=wildcards.sample_type,
            nucleotide=wildcards.nucleotide,
            taxid=wildcards.taxid
        )
    output: 
        fastq="{webinterface}/{start_date}_{run_id}_web/snakemake_results_{sample}/{sample_type}_{nucleotide}/organism_fastq/{taxid}.fastq.gz"
    shell:
        """
        ln -s {input.fastq} {output.fastq}
        """ 

rule link_filtered_fastq_kingdom:
    input: 
        fastq= lambda wildcards: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/kingdom_fastq/{kingdom}.fastq.gz".format(
            outdir=config['outdir'],
            start_date=wildcards.start_date,
            run_id=wildcards.run_id,
            sample=wildcards.sample,
            sample_type=wildcards.sample_type,
            nucleotide=wildcards.nucleotide,
            kingdom=wildcards.kingdom
        )
    output: 
        fastq="{webinterface}/{start_date}_{run_id}_web/snakemake_results_{sample}/{sample_type}_{nucleotide}/kingdom_fastq/{kingdom}.fastq.gz"
    shell:
        """
        ln -s {input.fastq} {output.fastq}
        """ 

rule link_filtered_fastq_unclassified:
    input: 
        fastq= lambda wildcards: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/unclassified_fastq/unclassified.fastq.gz".format(
            outdir=config['outdir'],
            start_date=wildcards.start_date,
            run_id=wildcards.run_id,
            sample=wildcards.sample,
            sample_type=wildcards.sample_type,
            nucleotide=wildcards.nucleotide
        )
    output: 
        fastq="{webinterface}/{start_date}_{run_id}_web/snakemake_results_{sample}/{sample_type}_{nucleotide}/unclassified_fastq/unclassified.fastq.gz"
    shell:
        """
        ln -s {input.fastq} {output.fastq}
        """ 


def filter_fastq_according_to_classification(wildcards):
    checkpoint_output_organism = checkpoints.tableview.get(**wildcards).output['organism_dir']
    organism_list = expand("{webinterface}/{start_date}_{run_id}_web/snakemake_results_{sample}/{sample_type}_{nucleotide}/organism_fastq/{taxid}.fastq.gz",
           webinterface=config['webinterface'],
           start_date=wildcards.start_date,
           run_id=wildcards.run_id,
           sample=wildcards.sample,
           sample_type=wildcards.sample_type,
           nucleotide=wildcards.nucleotide,
           taxid=glob_wildcards(os.path.join(checkpoint_output_organism, "organism_{taxid, \d+}.tsv")).taxid)

    checkpoint_output_kingdom = checkpoints.tableview.get(**wildcards).output['kingdom_dir']
    kingdom_list = expand("{webinterface}/{start_date}_{run_id}_web/snakemake_results_{sample}/{sample_type}_{nucleotide}/kingdom_fastq/{kingdom}.fastq.gz",
           webinterface=config['webinterface'],
           start_date=wildcards.start_date,
           run_id=wildcards.run_id,
           sample=wildcards.sample,
           sample_type=wildcards.sample_type,
           nucleotide=wildcards.nucleotide,
           kingdom=glob_wildcards(os.path.join(checkpoint_output_kingdom, "kingdom_{kingdom}.tsv")).kingdom)

    organism_kingdom_list = organism_list + kingdom_list
    return organism_kingdom_list

rule call_filter_fastqs:
    input:
        classified_fastqs=filter_fastq_according_to_classification,
        unclassified_fastqs=lambda wildcards: "{webinterface}/{start_date}_{run_id}_web/snakemake_results_{sample}/{sample_type}_{nucleotide}/unclassified_fastq/unclassified.fastq.gz".format(
            webinterface=config['webinterface'],
            start_date=wildcards.start_date,
            run_id=wildcards.run_id,
            sample=wildcards.sample,
            sample_type=wildcards.sample_type,
            nucleotide=wildcards.nucleotide
        )
    output:
        "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/fastq_filtering_done"
    shell:
        "touch {output}"

rule tableview_case:
    input: 
        case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/readcount_tableview.tsv".format(
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
        case="{outdir}/{start_date}_{run_id}/tableview/case_readcount_tableview.tsv"
    shell: 
        """
        cp {input.case} {output.case}
        """  

rule tableview_case_control:
    input: 
        case = lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/readcount_tableview.tsv".format(
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


rule detailed_stats_case:
    input: 
        raw_reads_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage1/samples/count_raw_reads.txt".format(
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
        trimmed_reads_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage1/trimming/count_bbduk_trimmed_reads.txt".format(
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
        kmer_input_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage2/kmer_input/count_kmer_input.txt".format(
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
        kraken_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage3/kraken/count_kraken_filtered_classified.txt".format(
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
        kaiju_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage3/kaiju/count_kaiju_filtered_classified.txt".format(
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
        species_genus_higher_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage4/count_species_genus_higher.txt".format(
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
        kmer_doublets_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage4/count_doublets.txt".format(
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
        kmer_singletons_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage4/count_singletons.txt".format(
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
        detected_and_missing_slices_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage5/count_detected_missing_slices.txt".format(
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
        subset_blast_reads_taxids_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage6/count_reads_taxid_SubsetBLAST.txt".format(
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
        kmer_subset_blast_reads_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage6/count_kmer_SubsetBLAST.txt".format(
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
        nt_blast_reads_taxids_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage7/count_reads_taxid_BLASTnt.txt".format(
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
        kmer_subset_blast_nt_blast_reads_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage7/count_kmer_SubsetBLAST_BLASTnt.txt".format(
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
        classified_reads_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/classified_reads_mincount.tsv".format(
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
        detailed_stats_out="{outdir}/{start_date}_{run_id}/tableview/case_detailed_stats.tsv"
    params: 
        case_control=False
    conda: "../../conda/R_env.yaml"
    script:  "../../scripts/reformat_results/detailed_stats.R"  

rule detailed_stats_case_control:
    input: 
        # Case
        raw_reads_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage1/samples/count_raw_reads.txt".format(
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
        trimmed_reads_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage1/trimming/count_bbduk_trimmed_reads.txt".format(
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
        kmer_input_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage2/kmer_input/count_kmer_input.txt".format(
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
        kraken_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage3/kraken/count_kraken_filtered_classified.txt".format(
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
        kaiju_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage3/kaiju/count_kaiju_filtered_classified.txt".format(
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
        species_genus_higher_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage4/count_species_genus_higher.txt".format(
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
        kmer_doublets_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage4/count_doublets.txt".format(
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
        kmer_singletons_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage4/count_singletons.txt".format(
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
        detected_and_missing_slices_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage5/count_detected_missing_slices.txt".format(
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
        subset_blast_reads_taxids_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage6/count_reads_taxid_SubsetBLAST.txt".format(
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
        kmer_subset_blast_reads_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage6/count_kmer_SubsetBLAST.txt".format(
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
        nt_blast_reads_taxids_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage7/count_reads_taxid_BLASTnt.txt".format(
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
        kmer_subset_blast_nt_blast_reads_case=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage7/count_kmer_SubsetBLAST_BLASTnt.txt".format(
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
        classified_reads_case = lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/classified_reads_mincount.tsv".format(
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

        # Control
        raw_reads_control=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage1/samples/count_raw_reads.txt".format(
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
        trimmed_reads_control=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage1/trimming/count_bbduk_trimmed_reads.txt".format(
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
        kmer_input_control=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage2/kmer_input/count_kmer_input.txt".format(
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
        kraken_control=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage3/kraken/count_kraken_filtered_classified.txt".format(
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
        kaiju_control=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage3/kaiju/count_kaiju_filtered_classified.txt".format(
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
        species_genus_higher_control=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage4/count_species_genus_higher.txt".format(
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
        kmer_doublets_control=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage4/count_doublets.txt".format(
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
        kmer_singletons_control=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage4/count_singletons.txt".format(
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
        detected_and_missing_slices_control=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage5/count_detected_missing_slices.txt".format(
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
        subset_blast_reads_taxids_control=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage6/count_reads_taxid_SubsetBLAST.txt".format(
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
        kmer_subset_blast_reads_control=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage6/count_kmer_SubsetBLAST.txt".format(
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
        nt_blast_reads_taxids_control=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage7/count_reads_taxid_BLASTnt.txt".format(
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
        kmer_subset_blast_nt_blast_reads_control=lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage7/count_kmer_SubsetBLAST_BLASTnt.txt".format(
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
        classified_reads_control = lambda wildcards: "{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/classified_reads_mincount.tsv".format(
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
        detailed_stats_out="{outdir}/{start_date}_{run_id}/tableview/case_control_detailed_stats.tsv"
    params: 
        case_control=True
    conda: "../../conda/R_env.yaml"
    script: "../../scripts/reformat_results/detailed_stats.R"   



