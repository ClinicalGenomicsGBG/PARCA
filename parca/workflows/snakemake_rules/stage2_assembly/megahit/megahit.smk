
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
        "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_RNA/stage1/fiona/trimmed_reads_fiona.fa"
    output:
        done_file=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_RNA/stage2/megahit/done"),
        intermediate_contigs_dir=temp(directory("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_RNA/stage2/megahit/intermediate_contigs")),
        opts=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_RNA/stage2/megahit/options.json"),
        contigs= temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_RNA/stage2/megahit/RNA.contigs.fa"),
        log=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_RNA/stage2/megahit/RNA.log"),
        read_count= temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_SE_RNA/stage2/megahit/count_assembled_contigs.txt")
    params:     
        out_prefix="RNA",
        outdir="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_RNA/stage2/megahit",
        min_contig_length=100
    # conda: "../../../conda/bbmap_env.yaml" #config['conda_environment']
    singularity: config['singularity_bbmap_env']
    log: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/logs_SE_RNA/stage2/megahit.log"
    benchmark: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/benchmarks_SE_RNA/stage2/megahit.log"
    threads: 110
    shell: 
        """
        megahit \
            -t {threads} \
            --out-dir {params.outdir}/tmp \
            --read {input} \
            --out-prefix {params.out_prefix} \
            --min-contig-len {params.min_contig_length} \
            &> {log};
        mv {params.outdir}/tmp/* {params.outdir};
        rmdir {params.outdir}/tmp;
        echo $(grep ">" {output.contigs}|wc -l) > {output.read_count};
        """
        #echo $(cat {output.contigs}|wc -l)/4|bc  > {output.read_count};


rule megahit_PE_RNA:
    """ 
    Rule for running Megahit metagenomic assembler on paired end RNA.
    Input: 
        unmerged=interleaved paired end reads.
        merged=merged paired end reads.
    Params: 
        out_prefix=Output prefix.
        outdir=Output directory.
        min_contig_len=The minimum length for being assigned as a contig.
    Output: 
        Log files, stats and the contigs file in Fasta format.
    """ 
    input:
        unmerged="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_RNA/stage1/trimming/unmerged_reads_trimmed.fq",
        merged="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_RNA/stage1/trimming/merged_reads_trimmed.fq"
    output:
        done_file=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_RNA/stage2/megahit/done"),
        intermediate_contigs_dir=temp(directory("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_RNA/stage2/megahit/intermediate_contigs")),
        opts=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_RNA/stage2/megahit/options.json"),
        contigs=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_RNA/stage2/megahit/RNA.contigs.fa"),
        log=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_RNA/stage2/megahit/RNA.log"),
        read_count=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_PE_RNA/stage2/megahit/count_assembled_contigs.txt")
    params:     
        out_prefix="RNA",
        outdir="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_RNA/stage2/megahit",
        #sub_outdir="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_RNA/stage2",
        min_contig_length=100
    # conda: "../../../conda/bbmap_env.yaml" #config['conda_environment']
    singularity: config['singularity_bbmap_env']
    log: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/logs_PE_RNA/stage2/megahit.log"
    benchmark: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/benchmarks_PE_RNA/stage2/megahit.log"
    threads: 110
    shell: 
        """
        megahit \
            -t {threads} \
            --out-dir {params.outdir}/tmp \
            --12 {input.unmerged} \
            --read {input.merged} \
            --out-prefix {params.out_prefix} \
            --min-contig-len {params.min_contig_length} \
            &> {log};
        mv {params.outdir}/tmp/* {params.outdir};
        rmdir {params.outdir}/tmp;
        echo $(grep ">" {output.contigs}|wc -l) > {output.read_count};
        """


# "megahit -t 110 --out-dir $outdir/megahit $inreads --out-prefix RNA $mincontig"
#-t/--num-cpu-threads     <int>          number of CPU threads, at least 2 if GPU enabled. [# of logical processors]
#-o/--out-dir             <string>       output directory [./megahit_out]
#-r/--read                <se>           comma-separated list of fasta/q single-end files
#--out-prefix             <string>       output prefix (the contig file will be OUT_DIR/OUT_PREFIX.contigs.fa)
#--min-contig-len         <int>          minimum length of contigs to output [200]
# --12                     <pe12>         comma-separated list of interleaved fasta/q paired-end files
# -r/--read                <se>           comma-separated list of fasta/q single-end files

# $assemblystatsline= "cat $outdir/megahit/RNA.contigs.fa|grep '>'|wc -l";
# $contigcount=`$assemblystatsline`;
# chomp $contigcount;
# $statisticshash{$sampletype}{"Assembled Contigs"}=$contigcount;