
rule readcount_RNA:
    input: 
        cov="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_RNA/stage2/pileup/bbmap_cov.txt"
    output: 
        readcount="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_RNA/stage8/readcount.tsv"
    script: "../../scripts/readcount_formatting.R" 

rule generate_krona_plot:
    input:
        ""
    output:
        ""
    conda: "../../conda/krona.yaml"
    shell:
        """
        ktImportText $kronaline -o $mainoutdir/text.krona.html
        """