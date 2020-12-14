
rule readcount_RNA:
    input: 
        cov="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_RNA/stage2/pileup/bbmap_cov.txt"
    output: 
    script: "../../scripts/readcount_formatting.R" 