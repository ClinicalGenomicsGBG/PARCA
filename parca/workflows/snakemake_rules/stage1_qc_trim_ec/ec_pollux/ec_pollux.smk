
rule pollux:
    input: 
        "{outdir}/snakemake_results_{sample}/stage1/trimming/trimmed_reads.fq"
    output: 
        "{outdir}/snakemake_results_{sample}/stage1/pollux/trimmed_reads.corrected.fq"
    params:
        pollux=config['pollux_path'],
        outdir="{outdir}/snakemake_results_{sample}/stage1/pollux"
        corrected="{outdir}/snakemake_results_{sample}/stage1/pollux/trimmed_reads.fq.corrected"
    shell:
        """
        {params.pollux} \ 
            -s false \ 
            -n true \ 
            -d true \ 
            -h true \ 
            -f false \ 
            -i {input} \ 
            -o {params.outdir}; \
        mv {params.corrected} {output};
        """

# system ("pollux -s false -n true -d true -h true -f false -i $outdir/trimmed_reads.fq -o $outdir/pollux")
# system("mv $outdir/pollux/trimmed_reads.fq.corrected $outdir/pollux/trimmed_reads.corrected.fq")