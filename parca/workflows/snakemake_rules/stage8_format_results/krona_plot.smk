
rule readcount_RNA:
    input: 
        cov="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_RNA/stage2/pileup/bbmap_cov.txt",
        all_classed_read_taxid_names="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_RNA/stage8/all_classed_read_taxid_names.txt"
    output: 
        readcount="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_RNA/stage8/readcount.tsv",
        readcount_krona="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_RNA/stage8/readcount_krona.tsv"
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
        readcount_krona="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_DNA/stage8/readcount_krona.tsv"
    params: 
        mincount=config['krona_plot_min_count'],
        DNA_or_RNA="DNA"
    conda: "../../conda/R_env.yaml"
    script: "../../scripts/reformat_results/readcount_formatting.R" 

rule tableview_SE:
    input: 
        read_count="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage8/readcount.tsv",
        trimmed_read_count="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_SE_{nucleotide}/stage1/trimming/count_bbduk_trimmed_reads.txt"
    output: 
        tableview="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage8/tableview/readcount_tableview.tsv",
        read_classifications="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage8/tableview/read_classifications.tsv"
        organism_dir=directory("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage8/tableview/organism_dir")
    params:
        SE_or_PE="SE",
        mincount=config['tableview_min_count']
    conda: "../../conda/R_env.yaml"
    script: "../../scripts/reformat_results/tableview_splitting.R"

rule tableview_PE:
    input: 
        read_count="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage8/readcount.tsv"
        trimmed_read_count_unmerged="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_PE_{nucleotide}/stage1/trimming/count_bbduk_unmerged_reads_trimmed_raw.txt",
        trimmed_read_count_merged="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_PE_{nucleotide}/stage1/trimming/count_bbduk_merged_reads_trimmed.txt"
    output: 
        tableview="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage8/tableview/readcount_tableview.tsv",
        read_classifications="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage8/tableview/read_classifications.tsv"
        organism_dir=directory("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage8/tableview/organism_dir")
    params:
        SE_or_PE="PE",
        mincount=config['tableview_min_count']
    conda: "../../conda/R_env.yaml"
    script: "../../scripts/reformat_results/tableview_splitting.R"

# rule generate_krona_plot_case:
#     input:
#         readcount_krona="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/krona/readcount_krona.tsv"
#     output:
#         krona_html="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/krona/text.krona.html"
#     conda: "../../conda/krona.yaml"
#     shell:
#         """
#         ktImportText {input.readcount_krona},"case" -o {output.krona_html};
#         """

# rule generate_krona_plot_case_control:
#     input:
#         readcount_krona_case="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/krona/readcount_krona.tsv"
#         readcount_krona_control="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/krona/readcount_krona.tsv"
#     output:
#         krona_html="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/krona/text.krona.html"
#     conda: "../../conda/krona.yaml"
#     shell:
#         """
#         ktImportText {input.readcount_krona_case},"case" {input.readcount_krona_control},"control" -o {output.krona_html};
#         """