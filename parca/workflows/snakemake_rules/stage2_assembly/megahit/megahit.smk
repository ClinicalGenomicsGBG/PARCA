
# done
# intermediate_contigs
# opts.txt
# RNA.contigs.fa
# RNA.log
rule megahit_SE_RNA:
    """ 
    Rule for running Megahit metagenomic assembler on single end RNA.
    Input: 
        A trimmed and error corrected fastq file.
    Params: 
        out_prefix=Output prefix.
        outdir=Output directory.
        min_contig_len=The minimum length for being assigned as a contig.
    Output: 
        Log files, stats and the contigs file in Fasta format.
    """ 
    input:
        "{outdir}/snakemake_results_{sample}/SE_RNA/stage1/fiona/trimmed_reads_fiona.fa"
    output:
        done_file="{outdir}/snakemake_results_{sample}/SE_RNA/stage2/megahit/done",
        intermediate_contigs_dir=directory("{outdir}/snakemake_results_{sample}/SE_RNA/stage2/megahit/intermediate_contigs"),
        opts="{outdir}/snakemake_results_{sample}/SE_RNA/stage2/megahit/options.json",
        contigs= "{outdir}/snakemake_results_{sample}/SE_RNA/stage2/megahit/RNA.contigs.fa",
        log="{outdir}/snakemake_results_{sample}/SE_RNA/stage2/megahit/RNA.log",
        read_count= "{outdir}/snakemake_results_{sample}/stats_SE_RNA/stage2/megahit/count_assembled_contigs.txt"
    params:     
        out_prefix="RNA",
        outdir="{outdir}/snakemake_results_{sample}/SE_RNA/stage2/megahit",
        min_contig_length=100
    conda: "../../../conda/bbmap_env.yaml" #config['conda_environment']
    log: "{outdir}/snakemake_results_{sample}/logs_SE_RNA/stage2/megahit.txt"
    benchmark: "{outdir}/snakemake_results_{sample}/benchmarks_SE_RNA/stage2/megahit.txt"
    threads: 110
    shell: 
        """
        megahit -t {threads} --out-dir {params.outdir} --read {input} --out-prefix {params.out_prefix} --min-contig-len {params.min_contig_length} &> {log};
        echo $(grep ">" {output.contigs}|wc -l) > {output.read_count};
        """
        #echo $(cat {output.contigs}|wc -l)/4|bc  > {output.read_count};

#-t/--num-cpu-threads     <int>          number of CPU threads, at least 2 if GPU enabled. [# of logical processors]
#-o/--out-dir             <string>       output directory [./megahit_out]
#-r/--read                <se>           comma-separated list of fasta/q single-end files
#--out-prefix             <string>       output prefix (the contig file will be OUT_DIR/OUT_PREFIX.contigs.fa)
#--min-contig-len         <int>          minimum length of contigs to output [200]

# $assemblystatsline= "cat $outdir/megahit/RNA.contigs.fa|grep '>'|wc -l";
# $contigcount=`$assemblystatsline`;
# chomp $contigcount;
# $statisticshash{$sampletype}{"Assembled Contigs"}=$contigcount;