
rule fiona_SE_RNA:
    input: 
        "{outdir}/snakemake_results_{sample}/SE_RNA/stage1/pollux/trimmed_reads.corrected.fq"
    output: 
        "{outdir}/snakemake_results_{sample}/SE_RNA/stage1/fiona/trimmed_reads_fiona.fq"
    params:
        fiona=config['fiona_path']
    log:  "{outdir}/snakemake_results_{sample}/logs_SE_RNA/stage1/fiona.log"
    benchmark: "{outdir}/snakemake_results_{sample}/benchmarks_SE_RNA/stage1/fiona.txt"
    threads: 110
    shell:
        """
        {params.fiona} -nt {threads} -g 100000000 {input} {output} &> {log}
        """

rule fiona_SE_DNA:
    input: 
       "{outdir}/snakemake_results_{sample}/SE_DNA/stage1/trimming/trimmed_reads.fq"
    output: 
        "{outdir}/snakemake_results_{sample}/SE_DNA/stage1/fiona/trimmed_reads_fiona.fq"
    params:
        fiona=config['fiona_path']
    log:  "{outdir}/snakemake_results_{sample}/logs_SE_DNA/stage1/fiona.log"
    benchmark: "{outdir}/snakemake_results_{sample}/benchmarks_SE_DNA/stage1/fiona.txt"
    threads: 110
    shell:
        """
        {params.fiona} -nt {threads} -g 100000000 {input} {output} &> {log}
        """

#fiona -nt 110 -g 100000000 $outdir/pollux/trimmed_reads.corrected.fq $outdir/trimmed_reads_fiona.fa

# -nt, --num-threads INTEGER
#     Number of threads to use (default 1). In range [1..inf]. Default: 1.
# -g, --genome-length INT64
#    Approximate length of the underlying genome. In range [1..inf].