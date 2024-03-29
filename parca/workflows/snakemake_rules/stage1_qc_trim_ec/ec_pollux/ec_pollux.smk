# Maintainer Pernilla Ericsson

rule pollux_SE_RNA:
    """ 
    Rule for error correction using Pollux.
    Input: 
        Trimmed fastq files.
    Params: 
        pollux=Path to pollux software.
        outdir=Out directory path.
        corrected=intermediate file.
    Output: 
        Error corrected reads.
    """ 
    input: 
        "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_RNA/stage1/trimming/trimmed_reads.fq"
    output: 
        temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_RNA/stage1/pollux/trimmed_reads.corrected.fq")
    params:
        pollux=config['pollux_path'], #config['pollux_path'],
        outdir="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_RNA/stage1/pollux",
        corrected="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_RNA/stage1/pollux/trimmed_reads.fq.corrected"
    log: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/logs_SE_RNA/stage1/pollux.log"
    benchmark: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/benchmarks_SE_RNA/stage1/pollux.log"
    shell:
        """
        {params.pollux} \
            -s false \
            -n true \
            -d true \
            -h true \
            -f false \
            -i {input} \
            -o {params.outdir} \
            &> {log}; \
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