
rule bbduk_trimming:
    input:
        lambda wildcards: "{sampledir}/{sample}.fastq".format(
                            sampledir=config['sampledir'], 
                            sample=wildcards.sample)
    output:
        "{outdir}/snakemake_results_{sample}/stage1/trimming"
    shell:
        """
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