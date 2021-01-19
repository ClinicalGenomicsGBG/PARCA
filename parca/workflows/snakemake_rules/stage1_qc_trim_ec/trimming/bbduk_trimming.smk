
rule bbduk_trimming_SE:
    """ 
    Rule for trimming fastq files.
    Input: 
        A fastq file.
    Params: 
        adapters=adapters to be removed.
        adaptertrimcommand=trimming settings.
    Output: 
        reads=Trimmed reads.
        stats=Statistics generated from bbduk.
        trimmed_read_count=The number of reads after trimming.
    """ 
    input:
        "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/samples/{sample}.fastq"
    output:
        reads= temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage1/trimming/trimmed_reads.fq"),
        stats= temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_SE_{nucleotide}/stage1/trimming/bbduk_stats.txt"),
        trimmed_read_count=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_SE_{nucleotide}/stage1/trimming/count_bbduk_trimmed_reads.txt")
    params:     
        adapters=','.join(list(metadata_df.loc[metadata_df['sample_id'] == '{sample}']['adapters'])),
        adaptertrimcommand="ktrim=l k=16 mink=11 hdist=1 rcomp=t"
    #threads: 23
    #conda: "../../../conda/bbmap_env.yaml" #config['conda_environment']
    singularity: config['singularity_bbmap_env']
    log: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/logs_SE_{nucleotide}/stage1/trimming.log"
    benchmark: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/benchmarks_SE_{nucleotide}/stage1/trimming.txt"
    shell:
        """
        adapter={params.adapters};
        if [[ -z ${{adapter//NA/}} ]] || [[ ${{adapter//NA/}} == "," ]]; then
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
            adapter=$(echo $adapter | sed 's/^,//g' | sed 's/,$//g');
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
                ref=$adapter \
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
        interleaved="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/samples/{sample}_interleaved.fastq"
    output: 
        merged=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/merged_reads.fastq"),
        unmerged=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/unmerged_reads.fastq")
    params:     
        mininsert=17
    # conda: "../../../conda/bbmap_env.yaml" #config['conda_environment']
    singularity: config['singularity_bbmap_env']
    log: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/logs_PE_{nucleotide}/stage1/bbduk_merging.log"
    benchmark: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/benchmarks_PE_{nucleotide}/stage1/bbduk_merging.txt"
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
    Rule for trimming merged fastq files with parameter lKtrim for adapters.
    Input: 
        A fastq file.
    Params: 
        adapters=adapters to be removed.
        adaptertrimcommand=trimming settings.
    Output: 
        reads=Trimmed reads.
        stats=Statistics generated from bbduk.
        trimmed_read_count=The number of reads after trimming.
    """ 
    input:
        "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/merged_reads.fastq"
    output:
        reads=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/merged_reads_lKtrim.fq"),
        stats=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_PE_{nucleotide}/stage1/trimming/bbduk_stats_merged_reads_lKtrim.txt"),
        trimmed_read_count=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_PE_{nucleotide}/stage1/trimming/count_bbduk_merged_reads_lKtrim.txt")
    params:     
        adapters=','.join(list(metadata_df.loc[metadata_df['sample_id'] == '{sample}']['adapters'])),
        adaptertrimcommand="ktrim=l k=16 mink=11 hdist=1 rcomp=t"
    # conda: "../../../conda/bbmap_env.yaml" #config['conda_environment']
    singularity: config['singularity_bbmap_env']
    log: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/logs_PE_{nucleotide}/stage1/merged_reads_lKtrim.log"
    benchmark: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/benchmarks_PE_{nucleotide}/stage1/merged_reads_lKtrim.txt"
    shell:
        """
        adapter={params.adapters};
        if [[ -z ${{adapter//NA/}} ]] || [[ ${{adapter//NA/}} == "," ]]; then
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
            adapter=$(echo $adapter | sed 's/^,//g' | sed 's/,$//g');
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
                ref=$adapter \
                {params.adaptertrimcommand} \
                &> {log}; \
        fi; \
        echo count > {output.trimmed_read_count};
        echo $(cat {output.reads}|wc -l)/4|bc  >> {output.trimmed_read_count};
        """


rule bbduk_trimming_PE_merged_rKtrim:
    """ 
    Rule for trimming merged fastq files with parameter rKtrim for adapters.
    Input: 
        A fastq file.
    Params: 
        adapters=adapter to be removed.
        adaptertrimcommand=trimming settings.
    Output: 
        reads=Trimmed reads.
        stats=Statistics generated from bbduk.
        trimmed_read_count=The number of reads after trimming.
    """ 
    input:
        "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/merged_reads_lKtrim.fq"
    output:
        reads=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/merged_reads_trimmed.fq"),
        stats=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_PE_{nucleotide}/stage1/trimming/bbduk_stats_merged_reads_trimmed.txt"),
        trimmed_read_count=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_PE_{nucleotide}/stage1/trimming/count_bbduk_merged_reads_trimmed.txt")
    params:     
        adapters=','.join(list(metadata_df.loc[metadata_df['sample_id'] == '{sample}']['adapters'])),
        adaptertrimcommand="ktrim=r k=16 mink=11 hdist=1 rcomp=t" 
    # conda: "../../../conda/bbmap_env.yaml" #config['conda_environment']
    singularity: config['singularity_bbmap_env']
    log: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/logs_PE_{nucleotide}/stage1/merged_reads_trimmed.log"
    benchmark: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/benchmarks_PE_{nucleotide}/stage1/merged_reads_trimmed.txt"
    shell:
        """
        adapter={params.adapters};
        if [[ -z ${{adapter//NA/}} ]] || [[ ${{adapter//NA/}} == "," ]]; then
            bbduk.sh \
                in={input} \
                stats={output.stats} \
                out={output.reads} \
                minlength=40 \
                overwrite=true \
                &> {log}; \
        else
            adapter=$(echo $adapter | sed 's/^,//g' | sed 's/,$//g');
            bbduk.sh \
                in={input} \
                stats={output.stats} \
                out={output.reads} \
                minlength=40 \
                overwrite=true \
                ref=$adapter \
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
        adapters=adapter to be removed.
        adaptertrimcommand=trimming settings.
    Output: 
        reads=Trimmed reads.
        stats=Statistics generated from bbduk.
        trimmed_read_count=The number of reads after trimming.
    """ 
    input:
        "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/unmerged_reads.fastq"
    output:
        reads=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/unmerged_reads_trimmed_raw.fq"),
        stats=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_PE_{nucleotide}/stage1/trimming/bbduk_stats_unmerged_reads_trimmed_raw.txt"),
        trimmed_read_count=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_PE_{nucleotide}/stage1/trimming/count_bbduk_unmerged_reads_trimmed_raw.txt")
    params:     
        adapters=','.join(list(metadata_df.loc[metadata_df['sample_id'] == '{sample}']['adapters'])),
        adaptertrimcommand="ktrim=l k=16 mink=11 hdist=1 rcomp=t" 
    # conda: "../../../conda/bbmap_env.yaml" #config['conda_environment']
    singularity: config['singularity_bbmap_env']
    log: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/logs_PE_{nucleotide}/stage1/unmerged_reads_trimmed_raw.log"
    benchmark: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/benchmarks_PE_{nucleotide}/stage1/unmerged_reads_trimmed_raw.txt"
    shell:
        """
        adapter={params.adapters};
        if [[ -z ${{adapter//NA/}} ]] || [[ ${{adapter//NA/}} == "," ]]; then
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
            adapter=$(echo $adapter | sed 's/^,//g' | sed 's/,$//g');
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
                ref=$adapter \
                {params.adaptertrimcommand} \
                &> {log}; \
        fi; \
        echo count > {output.trimmed_read_count};
        echo $(cat {output.reads}|wc -l)/4|bc  >> {output.trimmed_read_count};
        """

rule reformat_unmerged_PE:
    input: 
        reads="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/unmerged_reads_trimmed_raw.fq"
    output: 
        reads=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/unmerged_reads_trimmed.fq")
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

rule total_trimmed_reads_PE:
    input:
        trimmed_read_count_unmerged="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_PE_{nucleotide}/stage1/trimming/count_bbduk_unmerged_reads_trimmed_raw.txt",
        trimmed_read_count_merged="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_PE_{nucleotide}/stage1/trimming/count_bbduk_merged_reads_trimmed.txt" 
    output:
        trimmed_read_count=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_PE_{nucleotide}/stage1/trimming/count_bbduk_trimmed_reads.txt") 
    shell:
        """
        echo count > {output.trimmed_read_count};
        echo $(grep "^[0-9]" {input.trimmed_read_count_unmerged})+$(grep "^[0-9]" {input.trimmed_read_count_merged})|bc >> {output.trimmed_read_count};
        """ 

 
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


