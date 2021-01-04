from ...utils.process_runinfo_metadata import ProcessRuninfoMetadata


rule readcount_RNA:
    input: 
        cov="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_RNA/stage2/pileup/bbmap_cov.txt",
        all_classed_read_taxid_names="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_RNA/stage8/all_classed_read_taxid_names.txt"
    output: 
        readcount="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_RNA/stage8/readcount.tsv",
        readcount_krona="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_RNA/stage8/krona/readcount_krona.tsv"
    params: 
        mincount=config['krona_plot_min_count'],
        DNA_or_RNA="RNA"
    conda: "../../conda/R_env.yaml"
    script: "../../scripts/reformat_results/readcount_formatting.R" 

rule readcount_DNA:
    input: 
        all_classed_read_taxid_names="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_DNA/stage8/all_classed_read_taxid_names.txt",
    output: 
        readcount="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_DNA/stage8/readcount.tsv",
        readcount_krona="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_DNA/stage8/krona/readcount_krona.tsv"
    params: 
        mincount=config['krona_plot_min_count'],
        DNA_or_RNA="DNA"
    conda: "../../conda/R_env.yaml"
    script: "../../scripts/reformat_results/readcount_formatting.R" 


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