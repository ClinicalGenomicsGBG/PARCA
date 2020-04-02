
# done
# intermediate_contigs
# opts.txt
# RNA.contigs.fa
# RNA.log
rule megahit_SE_RNA:
    input:
        "{outdir}/snakemake_results_{sample}/SE_RNA/stage1/fiona/trimmed_reads_fiona.fq"
    output:
        done_file="{outdir}/snakemake_results_{sample}/SE_RNA/stage2/megahit/done",
        intermediate_contigs_dir=directory("{outdir}/snakemake_results_{sample}/SE_RNA/stage2/megahit/intermediate_contigs"),
        opts="{outdir}/snakemake_results_{sample}/SE_RNA/stage2/megahit/opts.txt",
        contigs= "{outdir}/snakemake_results_{sample}/SE_RNA/stage2/megahit/RNA.contigs.fa",
        log="{outdir}/snakemake_results_{sample}/SE_RNA/stage2/megahit/RNA.log",
        read_count= "{outdir}/snakemake_results_{sample}/stats_SE_RNA/stage2/megahit/count_assembled_contigs.txt"
    params:     
        out_prefix="RNA",
        outdir="{outdir}/snakemake_results_{sample}/SE_RNA/stage2/megahit"
    conda: config['conda_environment']
    log: "{outdir}/snakemake_results_{sample}/logs_SE_RNA/stage2/megahit.txt"
    benchmark: "{outdir}/snakemake_results_{sample}/benchmarks_SE_RNA/stage2/megahit.txt"
    threads: 110
    shell: 
        """
        megahit -f -t {threads} --out-dir {params.outdir} --read {input} --out-prefix {params.out_prefix} --min-contig-len 100 &> {log};
        echo $(cat {output.contigs}|wc -l)/4|bc  > {output.read_count};
        """

#-t/--num-cpu-threads     <int>          number of CPU threads, at least 2 if GPU enabled. [# of logical processors]
#-o/--out-dir             <string>       output directory [./megahit_out]
#-r/--read                <se>           comma-separated list of fasta/q single-end files
#--out-prefix             <string>       output prefix (the contig file will be OUT_DIR/OUT_PREFIX.contigs.fa)
#--min-contig-len         <int>          minimum length of contigs to output [200]

# $assemblystatsline= "cat $outdir/megahit/RNA.contigs.fa|grep '>'|wc -l";
# $contigcount=`$assemblystatsline`;
# chomp $contigcount;
# $statisticshash{$sampletype}{"Assembled Contigs"}=$contigcount;