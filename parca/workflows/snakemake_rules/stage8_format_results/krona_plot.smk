
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

rule name:
    input: 
    output: 
    script: "../../scripts/reformat_results/tableview_case_control" 


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