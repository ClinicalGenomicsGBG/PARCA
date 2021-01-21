# Maintainer Pernilla Ericsson

rule fiona_SE_RNA:
    """ 
    Rule for running fiona error correction on single end RNA.
    Input: 
        Trimmed reads.
    Params: 
        fiona=Path to fiona software.
    Output: 
        Error corrected reads.
    """ 
    input: 
        "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_RNA/stage1/pollux/trimmed_reads.corrected.fq"
    output: 
        temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_RNA/stage1/fiona/trimmed_reads_fiona.fa")
    params:
        fiona=config['fiona_path'] #config['fiona_path']
    log:  "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/logs_SE_RNA/stage1/fiona.log"
    benchmark: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/benchmarks_SE_RNA/stage1/fiona.log"
    threads: 110
    shell:
        """
        {params.fiona} \
            -nt {threads} \
            -g 100000000 \
            {input} \
            {output} &> {log};
        """

rule fiona_SE_DNA:
    """ 
    Rule for running fiona error correction on single end DNA.
    Input: 
        Trimmed reads.
    Params: 
        fiona=Path to fiona software.
    Output: 
        Error corrected reads.
    """ 
    input: 
       "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_DNA/stage1/trimming/trimmed_reads.fq"
    output: 
        temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_DNA/stage1/fiona/trimmed_reads_fiona.fa")
    params:
        fiona=config['fiona_path'] #config['fiona_path']
    log:  "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/logs_SE_DNA/stage1/fiona.log"
    benchmark: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/benchmarks_SE_DNA/stage1/fiona.log"
    threads: 110
    shell:
        """
        {params.fiona} \
            -nt {threads} \
            -id 3 \
            -g 4000000000 \
            {input} \
            {output} &> {log};
        """


#fiona -nt 110 -g 100000000 $outdir/pollux/trimmed_reads.corrected.fq $outdir/trimmed_reads_fiona.fa

# -nt, --num-threads INTEGER
#     Number of threads to use (default 1). In range [1..inf]. Default: 1.
# -g, --genome-length INT64
#    Approximate length of the underlying genome. In range [1..inf].
# -id, --indel-length INTEGER
#         Maximal indel length. Use 0 for correcting only substitutions and 1
#         for edit distance corrections on Illumina reads. In range [0..4].
#         Default: 1.