
rule bbduk_trimming_SE:
    """ 
    Rule for trimming fastq files.
    Input: 
        A fastq file.
    Params: 
        adapters=adapter-string if adapters should be removed.
        adaptertrimcommand=trimming settings.
    Output: 
        reads=Trimmed reads.
        stats=Statistics generated from bbduk.
        trimmed_read_count=The number of reads after trimming.
    """ 
    input:
        "{outdir}/snakemake_results_{sample}/SE_{nucleotide}/samples/{sample}.fastq"
    output:
        reads= "{outdir}/snakemake_results_{sample}/SE_{nucleotide}/stage1/trimming/trimmed_reads.fq",
        stats= "{outdir}/snakemake_results_{sample}/stats_SE_{nucleotide}/stage1/trimming/bbduk_stats.txt",
        trimmed_read_count="{outdir}/snakemake_results_{sample}/stats_SE_{nucleotide}/stage1/trimming/count_bbduk_trimmed_reads.txt"
    params:     
        adapters= runinfo_dict['adapters'], #config['adapters'],
        adaptertrimcommand="ktrim=l k=16 mink=11 hdist=1 rcomp=t"
    #threads: 23
    conda: "../../../conda/bbmap_env.yaml" #config['conda_environment']
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

rule bbduk_merging_PE:
    """ 
    Rule for merging PE reads based on overlap of mininsert size.
    Input: 
        interleaved=Interleaved file with PE reads.
    Params: 
        mininsert=Minimum insert size to merge reads.
    Output: 
        merged=Reads that could be merged.
        unmerged=Reads that were not merged.
    """ 
    input: 
        interleaved="{outdir}/snakemake_results_{sample}/PE_{nucleotide}/samples/{sample}_interleaved.fastq"
    output: 
        merged="{outdir}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/merged_reads.fastq",
        unmerged="{outdir}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/unmerged_reads.fastq"
    params:     
        mininsert=17
    conda: "../../../conda/bbmap_env.yaml" #config['conda_environment']
    log: "{outdir}/snakemake_results_{sample}/logs_PE_{nucleotide}/stage1/merging.log"
    benchmark: "{outdir}/snakemake_results_{sample}/benchmarks_PE_{nucleotide}/stage1/merging.txt"
    shell:
        """
        bbmerge.sh \
            mininsert={params.mininsert} \
            in={input.interleaved} \
            out={output.merged} \
            outu={output.unmerged} \
            &> {log};
        """ 

rule bbduk_trimming_PE_merged_lKtrim:
    """ 
    Rule for trimming fastq files with parameter lKtrim for adapters.
    Input: 
        A fastq file.
    Params: 
        adapters=adapter-string if adapters should be removed.
        adaptertrimcommand=trimming settings.
    Output: 
        reads=Trimmed reads.
        stats=Statistics generated from bbduk.
        trimmed_read_count=The number of reads after trimming.
    """ 
    input:
        "{outdir}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/merged_reads.fastq"
    output:
        reads="{outdir}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/merged_reads_lKtrim.fq",
        stats="{outdir}/snakemake_results_{sample}/stats_PE_{nucleotide}/stage1/trimming/bbduk_stats_merged_lKtrim.txt",
        trimmed_read_count="{outdir}/snakemake_results_{sample}/stats_PE_{nucleotide}/stage1/trimming/count_bbduk_merged_reads_lKtrim.txt"
    params:     
        adapters=runinfo_dict['adapters'], #config['adapters'],
        adaptertrimcommand="ktrim=l k=16 mink=11 hdist=1 rcomp=t"
    conda: "../../../conda/bbmap_env.yaml" #config['conda_environment']
    log: "{outdir}/snakemake_results_{sample}/logs_PE_{nucleotide}/stage1/trimming_merged_lKtrim.log"
    benchmark: "{outdir}/snakemake_results_{sample}/benchmarks_PE_{nucleotide}/stage1/trimming_merged_lKtrim.txt"
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
                trimq=12 \
                minlength=40 \
                qtrim=rl \
                overwrite=true \
                &> {log}; \
        else
            bbduk.sh \
                in={input} \
                entropymask=t \
                stats={output.stats} \
                out={output.reads} \
                entropy=0.9 \
                trimq=12 \
                minlength=40 \
                qtrim=rl \
                overwrite=true \
                ref={params.adapters} \
                {params.adaptertrimcommand} \
                &> {log}; \
        fi; \
        echo count > {output.trimmed_read_count};
        echo $(cat {output.reads}|wc -l)/4|bc  >> {output.trimmed_read_count};
        """


rule bbduk_trimming_PE_merged_rKtrim:
    """ 
    Rule for trimming fastq files with parameter rKtrim for adapters.
    Input: 
        A fastq file.
    Params: 
        adapters=adapter-string if adapters should be removed.
        adaptertrimcommand=trimming settings.
    Output: 
        reads=Trimmed reads.
        stats=Statistics generated from bbduk.
        trimmed_read_count=The number of reads after trimming.
    """ 
    input:
        "{outdir}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/trimmed_reads_lKtrim.fq"
    output:
        reads="{outdir}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/trimmed_reads_merged.fq",
        stats="{outdir}/snakemake_results_{sample}/stats_PE_{nucleotide}/stage1/trimming/bbduk_stats_merged.txt",
        trimmed_read_count="{outdir}/snakemake_results_{sample}/stats_PE_{nucleotide}/stage1/trimming/count_bbduk_trimmed_reads_merged.txt"
    params:     
        adapters=runinfo_dict['adapters'], #config['adapters'],
        adaptertrimcommand="ktrim=r k=16 mink=11 hdist=1 rcomp=t" 
    conda: "../../../conda/bbmap_env.yaml" #config['conda_environment']
    log: "{outdir}/snakemake_results_{sample}/logs_PE_{nucleotide}/stage1/trimming_merged.log"
    benchmark: "{outdir}/snakemake_results_{sample}/benchmarks_PE_{nucleotide}/stage1/trimming_merged.txt"
    shell:
        """
        adapter={params.adapters};
        if [[ -z ${{adapter/None/}} ]]; then
            bbduk.sh \
                in={input} \
                stats={output.stats} \
                out={output.reads} \
                minlength=40 \
                overwrite=true \
                &> {log}; \
        else
            bbduk.sh \
                in={input} \
                stats={output.stats} \
                out={output.reads} \
                minlength=40 \
                overwrite=true \
                ref={params.adapters} \
                {params.adaptertrimcommand} \
                &> {log}; \
        fi; \
        echo count > {output.trimmed_read_count};
        echo $(cat {output.reads}|wc -l)/4|bc  >> {output.trimmed_read_count};
        """

rule bbduk_trimming_PE_unmerged:
    """ 
    Rule for trimming unmerged paired fastq files.
    Input: 
        A fastq file.
    Params: 
        adapters=adapter-string if adapters should be removed.
        adaptertrimcommand=trimming settings.
    Output: 
        reads=Trimmed reads.
        stats=Statistics generated from bbduk.
        trimmed_read_count=The number of reads after trimming.
    """ 
    input:
        "{outdir}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/paired_reads.fastq"
    output:
        reads="{outdir}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/trimmed_reads_paired.fq",
        stats="{outdir}/snakemake_results_{sample}/stats_PE_{nucleotide}/stage1/trimming/bbduk_stats_paired.txt",
        trimmed_read_count="{outdir}/snakemake_results_{sample}/stats_PE_{nucleotide}/stage1/trimming/count_bbduk_trimmed_reads_paired.txt"
    params:     
        adapters=runinfo_dict['adapters'], #config['adapters'],
        adaptertrimcommand="ktrim=l k=16 mink=11 hdist=1 rcomp=t" 
    conda: "../../../conda/bbmap_env.yaml" #config['conda_environment']
    log: "{outdir}/snakemake_results_{sample}/logs_PE_{nucleotide}/stage1/trimming_paired.log"
    benchmark: "{outdir}/snakemake_results_{sample}/benchmarks_PE_{nucleotide}/stage1/trimming_paired.txt"
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
                trimq=12 \
                minlength=25 \
                qtrim=rl \
                overwrite=true \
                removeifeitherbad=f \
                &> {log}; \
        else
            bbduk.sh \
                in={input} \
                entropymask=t \
                stats={output.stats} \
                out={output.reads} \
                entropy=0.9 \
                trimq=12 \
                minlength=25 \
                qtrim=rl \
                overwrite=true \
                removeifeitherbad=f \
                ref={params.adapters} \
                {params.adaptertrimcommand} \
                &> {log}; \
        fi; \
        echo count > {output.trimmed_read_count};
        echo $(cat {output.reads}|wc -l)/4|bc  >> {output.trimmed_read_count};
        """

rule reformat_unmerged_PE:
    input: 
        reads="{outdir}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/trimmed_reads_paired.fq"
    output: 
        reads="{outdir}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/trimmed_reads_paired_reformatted.fq"
    run: 
        reformatted=[]
        with open(input['reads'], 'r') as filehandle:
            contents = filehandle.readlines()
            
            for i in range(0,len(contents),4):
                header=contents[i]
                sequence=contents[i+1]
                plus=contents[i+2]
                qual=contents[i+3]
                
                if header == "\n":
                    header="N\n"
                    qual="A\n"
                
                reformatted.append(header) 
                reformatted.append(sequence)
                reformatted.append(plus)
                reformatted.append(qual) 

        with open(output['reads'], 'w') as printresults:
            printresults.writelines("%s" % line for line in reformatted)

# rule count_trimmed_reads_PE:
#     input: 
#         "{outdir}/snakemake_results_{sample}/stats_PE_{nucleotide}/stage1/trimming/count_bbduk_trimmed_reads_paired.txt",
#         "{outdir}/snakemake_results_{sample}/stats_PE_{nucleotide}/stage1/trimming/count_bbduk_trimmed_reads_merged.txt"
#     output: 
#     shell:
#         """
#         """

 
# bbduk.sh
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

# bbmerge.sh
# Description:  Merges paired reads into single reads by overlap detection.
# With sufficient coverage, can also merge nonoverlapping reads by kmer extension.
# Kmer modes requires much more memory, and should be used with the bbmerge-auto.sh script.
# Please read bbmap/docs/guides/BBMergeGuide.txt for more information.
# Usage for interleaved files:	bbmerge.sh in=<reads> out=<merged reads> outu=<unmerged reads>
# Usage for paired files:     	bbmerge.sh in1=<read1> in2=<read2> out=<merged reads> outu1=<unmerged1> outu2=<unmerged2>
# mininsert: Minimum insert size to merge reads.
# in: Primary input. 'in2' will specify a second file.
# out: File for merged reads. 'out2' will specify a second file.
# outu: File for unmerged reads. 'outu2' will specify a second file.

# Parca v4.3b params
# $adaptertrimcommand="ref=$adapterpath ktrim=l k=16 mink=11 hdist=1 rcomp=t";
# $adaptertrimcommand_merged="ref=$adapterpath ktrim=l k=16 mink=11 hdist=1 rcomp=t";
# $adaptertrimcommand_merged_rev="ref=$adapterpath ktrim=r k=16 mink=11 hdist=1 rcomp=t";
# system ("bbduk.sh in=$outdir/rawreads.fastq entropymask=t stats=$outdir/bbduk_stats.txt out=$outdir/trimmed_reads.fq entropy=0.9 trimq=16 minlength=40 qtrim=rl overwrite=true $adaptertrimcommand_merged");


