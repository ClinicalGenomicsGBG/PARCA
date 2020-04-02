
rule pollux_SE_RNA:
    input: 
        "{outdir}/snakemake_results_{sample}/SE_RNA/stage1/trimming/trimmed_reads.fq"
    output: 
        "{outdir}/snakemake_results_{sample}/SE_RNA/stage1/pollux/trimmed_reads.corrected.fq"
    params:
        pollux=config['pollux_path'],
        outdir="{outdir}/snakemake_results_{sample}/SE_RNA/stage1/pollux",
        corrected="{outdir}/snakemake_results_{sample}/SE_RNA/stage1/pollux/trimmed_reads.fq.corrected"
    log: "{outdir}/snakemake_results_{sample}/logs_SE_RNA/stage1/pollux.log"
    benchmark: "{outdir}/snakemake_results_{sample}/benchmarks_SE_RNA/stage1/pollux.txt"
    shell:
        """
        {params.pollux} -s false -n true -d true -h true -f false -i {input} -o {params.outdir} &> {log}; \
        mv {params.corrected} {output};
        """

# system ("pollux -s false -n true -d true -h true -f false -i $outdir/trimmed_reads.fq -o $outdir/pollux")
# system("mv $outdir/pollux/trimmed_reads.fq.corrected $outdir/pollux/trimmed_reads.corrected.fq")

#-s      [bool]  Substitution corrections. "true" or "false".
#-n      [bool]  Insertion corrections. "true" or "false".
#-d      [bool]  Deletion corrections. "true" or "false".
#-h      [bool]  Homopolymer corrections. "true" or "false".
#-f      [bool]  Low k-mer read filtering. "true" or "false".
#-i      [file]  Specify one or many FASTQ input files.
#-o              Output directory.