
rule bbduk_trimming_SE:
    input:
         "{outdir}/snakemake_results_{sample}/SE_{nucleotide}/samples/{sample}.fastq"
    output:
        reads= "{outdir}/snakemake_results_{sample}/SE_{nucleotide}/stage1/trimming/trimmed_reads.fq",
        stats= "{outdir}/snakemake_results_{sample}/stats_SE_{nucleotide}/stage1/trimming/bbduk_stats.txt",
        trimmed_read_count="{outdir}/snakemake_results_{sample}/stats_SE_{nucleotide}/stage1/trimming/count_bbduk_trimmed_reads.txt"
    params:     
        adapters=config['adapters'],
        adaptertrimcommand="ktrim=l k=16 mink=11 hdist=1 rcomp=t"
    #threads: 23
    conda: config['conda_environment']
    log: "{outdir}/snakemake_results_{sample}/logs_SE_{nucleotide}/stage1/trimming.log"
    benchmark: "{outdir}/snakemake_results_{sample}/benchmarks_SE_{nucleotide}/stage1/trimming.txt"
    shell:
        """
        adapter={params.adapters};
        if [[ -z ${{adapter/None/}} ]]; then
            bbduk.sh \
                in={input} \
                entropymask=t \
                stats={output.stats} \
                out={output.reads} \
                entropy=0.9 \
                trimq=16 \
                minlength=40 \
                qtrim=rl \
                overwrite=true &> {log}; \
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
                {params.adaptertrimcommand} &> {log}; \
        fi; \
        echo count > {output.trimmed_read_count};
        echo $(cat {output.reads}|wc -l)/4|bc  >> {output.trimmed_read_count};
        """

# in=<file>           Main input. in=stdin.fq will pipe from stdin.
# entropymask=f       Values:
#                       f:  Discard low-entropy sequences.
#                       t:  Mask low-entropy parts of sequences with N.
#                       lc: Change low-entropy parts of sequences to lowercase.
# stats=<file>        Write statistics about which contamininants were detected.
# out=<file>          (outnonmatch) Write reads here that do not contain
#                     kmers matching the database.  'out=stdout.fq' will pipe
#                     to standard out.
# entropy=-1          Set between 0 and 1 to filter reads with entropy below
#                     that value.  Higher is more stringent.
# trimq=6             Regions with average quality BELOW this will be trimmed,
#                     if qtrim is set to something other than f.
# minlength=10        (ml) Reads shorter than this after trimming will be
#                     discarded.  Pairs will be discarded if both are shorter.
# qtrim=f             Trim read ends to remove bases with quality below trimq.
#                     Performed AFTER looking for kmers.
#                     Values:
#                             rl (trim both ends),
#                             f (neither end),
#                             r (right end only),
#                             l (left end only),
#                             w (sliding window).
# overwrite=t         (ow) Grant permission to overwrite files.
# ref=<file,file>     Comma-delimited list of reference files.
#                     You can also use ref=phix, ref=adapters, or ref=artifacts.
# ktrim=f             Trim reads to remove bases matching reference kmers.
#                     Values:
#                             f (don't trim),
#                             r (trim to the right),
#                             l (trim to the left)
# k=27                Kmer length used for finding contaminants.  Contaminants
#                     shorter than k will not be found.  k must be at least 1.
# mink=0              Look for shorter kmers at read tips down to this length,
#                     when k-trimming or masking.  0 means disabled.  Enabling
#                     this will disable maskmiddle.
# hammingdistance=0   (hdist) Maximum Hamming distance for ref kmers (subs only).
#                     Memory use is proportional to (3*K)^hdist.
# rcomp=t             Look for reverse-complements of kmers in addition to
#                     forward kmers.


# $adaptertrimcommand="ref=$adapterpath ktrim=l k=16 mink=11 hdist=1 rcomp=t";
# $adaptertrimcommand_merged="ref=$adapterpath ktrim=l k=16 mink=11 hdist=1 rcomp=t";
# $adaptertrimcommand_merged_rev="ref=$adapterpath ktrim=r k=16 mink=11 hdist=1 rcomp=t";


