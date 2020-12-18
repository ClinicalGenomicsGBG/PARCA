
rule readcount_RNA:
    input: 
        cov="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_RNA/stage2/pileup/bbmap_cov.txt",
        all_classed_read_taxid_names="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_RNA/stage8/all_classed_read_taxid_names.txt"
    output: 
        readcount="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_RNA/stage8/readcount.tsv",
        readcount_krona="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_RNA/stage8/readcount_krona.tsv"
    params: 
        mincount=config['krona_plot_min_count']
    script: "../../scripts/readcount_formatting_RNA.R" 

rule readcount_DNA:
    input: 
        all_classed_read_taxid_names="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_DNA/stage8/all_classed_read_taxid_names.txt"
    output: 
        readcount="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_DNA/stage8/readcount.tsv",
        readcount_krona="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_DNA/stage8/readcount_krona.tsv"
    params: 
        mincount=config['krona_plot_min_count']
    script: "../../scripts/readcount_formatting_DNA.R" 

rule generate_krona_plot:
    input:
        readcount="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/krona/readcount_krona.tsv"
    output:
        krona_html="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/krona/text.krona.html"
    conda: "../../conda/krona.yaml"
    # params: 
    #     run_type="case"
    shell:
        """
        ktImportText {input.readcount} -o {output.krona_html};
        """
# rule percent_calculation:
#     input: 
#         # all count
#         read_count="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_SE_RNA/stage2/kmer_input/count_kmer_input.txt"
#     output: 
#     run: 