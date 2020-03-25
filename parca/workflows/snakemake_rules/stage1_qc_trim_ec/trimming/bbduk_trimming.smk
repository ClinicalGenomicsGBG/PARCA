
rule bbduk_trimming:
    input:
        "{outdir}/snakemake_results_{sample}/{sample}.fastq"
    output:
        reads="{outdir}/snakemake_results_{sample}/stage1/trimming/trimmed_reads.fq",
        stats="{outdir}/snakemake_results_{sample}/stats/stage1/trimming/bbduk_stats.txt",
        trimmed_read_count="{outdir}/snakemake_results_{sample}/stats/stage1/trimming/bbduk_trimmed_read_count.txt"
    params:     
        adapters=config['adapters'],
        adaptertrimcommand="ktrim=l k=16 mink=11 hdist=1 rcomp=t"
    #threads: 23
    conda: config['conda_environment']
    benchmark: "{outdir}/snakemake_results_{sample}/benchmarks/stage1/trimming.txt"
    shell:
        """
        if [[ {params.adapters} == "None" ]]; then
            bbduk.sh \
                in={input} \
                entropymask=t \
                stats={output.stats} \
                out={output.reads} \
                entropy=0.9 \
                trimq=16 \
                minlength=40 \
                qtrim=rl \
                overwrite=true; \
        else
            bbduk.sh \
                in={input} \
                entropymask=t \
                stats={output.stats} \
                out={output.reads} \
                entropy=0.9 \
                trimq=16 \
                minlength=40 \
                qtrim=rl \
                overwrite=true \
                ref={params.adapters} \
                {params.adaptertrimcommand}; \
        fi; \
        echo $(cat {output.reads}|wc -l)/4|bc  > {output.trimmed_read_count}
        """

# system ("bbduk.sh 
        #in=$outdir/rawreads.fastq 
        #entropymask=t 
        #stats=$outdir/bbduk_stats.txt 
        #out=$outdir/trimmed_reads.fq 
        #entropy=0.9 
        #trimq=16 
        #minlength=40 
        #qtrim=rl 
        #overwrite=true 
        #$adaptertrimcommand_merged");
# $trimmedreads="trimmed_reads.fq";
# @trimmedarray=("trimmed_reads.fq");

# $adaptertrimcommand="ref=$adapterpath ktrim=l k=16 mink=11 hdist=1 rcomp=t";
# $adaptertrimcommand_merged="ref=$adapterpath ktrim=l k=16 mink=11 hdist=1 rcomp=t";
# $adaptertrimcommand_merged_rev="ref=$adapterpath ktrim=r k=16 mink=11 hdist=1 rcomp=t";


