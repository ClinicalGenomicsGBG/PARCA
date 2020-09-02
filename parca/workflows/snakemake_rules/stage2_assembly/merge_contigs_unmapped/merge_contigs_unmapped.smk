rule create_kmer_classifier_input_SE_RNA:
    """ 
    Rule for merging fasta files with contigs and fasta files with unmapped reads.
    Input: 
        contigs=contigs generated from megahit.
        unmapped_reads=Reads that could not be mapped back to the contigs.
    Output: 
        kmer_input=Merged contigs and unmapped reads.
        read_count=The number of sequences in kmer_input.
    """ 
    input:
        contigs="{outdir}/snakemake_results_{sample}/SE_RNA/stage2/megahit/RNA.contigs.fa",
        unmapped_reads="{outdir}/snakemake_results_{sample}/SE_RNA/stage2/bbwrap_alignment/unmapped_reads.fasta"
    output:
        kmer_input="{outdir}/snakemake_results_{sample}/SE_RNA/stage2/kmer_input/kmer_input.fasta",
        read_count="{outdir}/snakemake_results_{sample}/stats_SE_RNA/stage2/kmer_input/count_kmer_input.txt"
    #conda: config['conda_environment']
    shell:
        """
        cat {input.contigs} {input.unmapped_reads} > {output.kmer_input};
        
        echo count > {output.read_count};
        echo $(grep ">" {output.kmer_input}|wc -l) >> {output.read_count};
        """
        #echo $(cat {output.kmer_input}|wc -l)/4|bc  >> {output.read_count};

rule create_kmer_classifier_input_SE_DNA:
    """ 
    Rule for renaming the preprocessed (bbduk_trimming and fiona error correction) file.
    Input: 
        fasta=Sequences that were trimmed using bbduk and corrected using fiona.
    Output: 
        kmer_input=All sequences.
        read_count=The number of sequences in kmer_input.
    """ 
    input:
        fasta="{outdir}/snakemake_results_{sample}/SE_DNA/stage1/fiona/trimmed_reads_fiona.fa"
    output:
        kmer_input="{outdir}/snakemake_results_{sample}/SE_RNA/stage2/kmer_input/kmer_input.fasta",
        read_count="{outdir}/snakemake_results_{sample}/stats_SE_RNA/stage2/kmer_input/count_kmer_input.txt"
    #conda: config['conda_environment']
    shell:
        """
        cp {input.fasta} {output.kmer_input};
        
        echo count > {output.read_count};
        echo $(grep ">" {output.kmer_input}|wc -l) >> {output.read_count};
        """

rule join_unmerged_PE:
    """ 
    Rule for joining paired end reads using "N" as separator.
    Input: fasta=Reads that were not merged with bbmerge nor mapped to contigs from megahit.
    Output: fasta=Paired reads joined using "N" as separator.
    """ 
    input:
        fasta="{outdir}/snakemake_results_{sample}/PE_RNA/stage2/bbwrap_alignment/unmerged_reads_unmapped.fasta"
    output:
        fasta="{outdir}/snakemake_results_{sample}/PE_RNA/stage2/bbwrap_alignment/unmerged_reads_unmapped_joined.fasta"
    run:
        reformatted=[]
        with open(input['fasta'], 'r') as filehandle:
            contents = [line.strip() for line in filehandle.readlines()]

            for i in range(0,len(contents),4):
                headerFwd=contents[i]
                sequenceFwd=contents[i+1]

                headerRev=contents[i+2]
                sequenceRev=contents[i+3]

                readIdFwd=headerFwd.split(" ")[0]
                readIdRev=headerRev.split(" ")[0]
                if readIdFwd==readIdRev:
                    reformatted.append(readIdFwd)
                    reformatted.append(sequenceFwd+"N"+sequenceRev)

        with open(output['fasta'], 'w') as printresults:
            printresults.writelines("%s\n" % line for line in reformatted)            

rule create_kmer_classifier_input_PE_RNA:
    input: 
        unmerged_reads_unmapped="{outdir}/snakemake_results_{sample}/PE_RNA/stage2/bbwrap_alignment/unmerged_reads_unmapped_joined.fasta",
        contigs="{outdir}/snakemake_results_{sample}/PE_RNA/stage2/megahit/RNA.contigs.fa",
        merged_reads_mapped="{outdir}/snakemake_results_{sample}/PE_RNA/stage2/bbwrap_alignment/merged_reads_unmapped.fasta"
    output: 
        kmer_input="{outdir}/snakemake_results_{sample}/PE_RNA/stage2/kmer_input/kmer_input.fasta",
        read_count="{outdir}/snakemake_results_{sample}/stats_PE_RNA/stage2/kmer_input/count_kmer_input.txt"
    shell:
        """
        cat {input} > {output.kmer_input};

        echo count > {output.read_count};
        echo $(grep ">" {output.kmer_input}|wc -l) >> {output.read_count};
        """ 