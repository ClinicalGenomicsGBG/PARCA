rule bbwrap_alignment_SE_RNA:
    input:
        ref="{outdir}/snakemake_results_{sample}/SE_RNA/stage2/megahit/RNA.contigs.fa",
        reads="{outdir}/snakemake_results_{sample}/SE_RNA/stage1/fiona/trimmed_reads_fiona.fq"
    output:
        mapped="{outdir}/snakemake_results_{sample}/SE_RNA/stage2/bbwrap_alignment/aln.sam.gz",
        unmapped="{outdir}/snakemake_results_{sample}/SE_RNA/stage2/bbwrap_alignment/unmapped_reads.fasta",
        scafstats="{outdir}/snakemake_results_{sample}/stats_SE_RNA/stage2/bbwrap_alignment/bbmap_scafstats.txt",
        stats="{outdir}/snakemake_results_{sample}/stats_SE_RNA/stage2/bbwrap_alignment/bbmap_stats.txt"
    conda: config['conda_environment']
    log: "{outdir}/snakemake_results_{sample}/logs_SE_RNA/stage2/bbwrap_alignment.log"
    benchmark: "{outdir}/snakemake_results_{sample}/benchmarks_SE_RNA/stage2/bbwrap.txt"
    threads: 110
    shell: 
        """
        bbwrap.sh \
            ref={input.ref} \
            nodisk \
            in={input.reads} \
            out={output.mapped} \
            kfilter=22 \
            32bit=t \
            threads={threads} \
            maxindel=80 \
            perfectmode=t \
            secondarycov=f \
            secondary=f \
            scafstats={output.scafstats} \
            statsfile={output.stats} \
            outu={output.unmapped} \
            fastawrap=10000 &> {log};
        """

rule pileup_SE_RNA:
    input:
        "{outdir}/snakemake_results_{sample}/SE_RNA/stage2/bbwrap_alignment/aln.sam.gz"
    output:
        "{outdir}/snakemake_results_{sample}/SE_RNA/stage2/pileup/bbmap_cov.txt"
    conda: config['conda_environment']
    log: "{outdir}/snakemake_results_{sample}/logs_SE_RNA/stage2/pileup.log"
    shell:
        """
        pileup.sh \
            32bit=t \
            in={input} \
            out={output} &> {log};
        """


# "bbwrap.sh ref=$outdir/megahit/RNA.contigs.fa nodisk in=$trimmedreads out=$outdir/aln.sam.gz kfilter=22 32bit=t threads=110 maxindel=80 perfectmode=t secondarycov=f secondary=f scafstats=$outdir/bbmap_stats.txt statsfile=$outdir/bbmap_stats2.txt outu=$outdir/unmapped_reads.fasta fastawrap=10000"
# ref=<file>              Specify the reference sequence.  Only do this ONCE,
#                         when building the index (unless using 'nodisk').
# nodisk=f                Set to true to build index in memory and write nothing
#                         to disk except output.
#  in=                    stdin will accept reads from standard in, and out=stdout will write to
#                         standard out, but file extensions are still needed to specify the format of the
#                         input and output files e.g. in=stdin.fa.gz will read gzipped fasta from
#                         standard in; out=stdout.sam.gz will write gzipped sam.
# out=<file>              Write all reads to this file.
# kfilter=0               If positive, potential mapping sites must have at
#                         least this many consecutive exact matches.
# 32bit=f                 Set to true if you need per-base coverage over 64k.
# maxindel=16000          Don't look for indels longer than this. Lower is faster.
#                         Set to >=100k for RNAseq with long introns like mammals.
# perfectmode=f           Allow only perfect mappings when set to true (very fast).
# secondarycov=t          Include coverage of secondary alignments.
# secondary=f             Print secondary alignments.
# scafstats=<file>        Statistics on how many reads mapped to which scaffold.
# statsfile=stderr        Mapping statistics are printed here.
# outu=<file>             Write only unmapped reads to this file.  Does not
#                         include unmapped paired reads with a mapped mate.
# fastawrap=10000;    

# "pileup.sh 32bit=t in=$outdir/aln.sam.gz out=$outdir/bbmap_cov.txt"
# 32bit=f             Set to true if you need per-base coverage over 64k; does not affect per-scaffold coverage precision.
# in=<file>           The input sam file; this is the only required parameter.
# out=<file>          (covstats) Per-scaffold coverage info.

#bbwrap.sh ref=/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage2/megahit/RNA.contigs.fa in=/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/samples/a.fastq out=mapped.sam append