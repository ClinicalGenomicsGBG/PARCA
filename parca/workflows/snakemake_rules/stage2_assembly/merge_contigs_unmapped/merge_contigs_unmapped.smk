import re

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
        contigs="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_RNA/stage2/megahit/RNA.contigs.fa",
        unmapped_reads="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_RNA/stage2/bbwrap_alignment/unmapped_reads.fasta"
    output:
        kmer_input=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_RNA/stage2/kmer_input/kmer_input.fasta"),
        read_count=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_SE_RNA/stage2/kmer_input/count_kmer_input.txt")
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
        fasta="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_DNA/stage1/fiona/trimmed_reads_fiona.fa"
    output:
        kmer_input=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_DNA/stage2/kmer_input/kmer_input.fasta"),
        read_count=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_SE_DNA/stage2/kmer_input/count_kmer_input.txt")
    shell:
        """
        cp {input.fasta} {output.kmer_input};
        
        echo count > {output.read_count};
        echo $(grep ">" {output.kmer_input}|wc -l) >> {output.read_count};
        """

rule join_unmerged_PE_RNA:
    """ 
    Rule for joining paired end reads using "N" as separator.
    Input: fasta=Reads that were not merged with bbmerge nor mapped to contigs from megahit.
    Output: fasta=Paired reads joined using "N" as separator.
    """ 
    input:
        fasta="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_RNA/stage2/bbwrap_alignment/unmerged_reads_unmapped.fasta"
    output:
        fasta=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_RNA/stage2/bbwrap_alignment/unmerged_reads_unmapped_joined.fasta")
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
    """ 
    Rule for merging the contigs generated from megahit and the reads that could not map to the contigs unmerged_reads_unmapped and merged_reads_unmapped.
    Input: 
        unmerged_reads_unmapped=Reads that were not merged by bbduk and that could not map to the contigs.
        contigs=The contigs created from megahit.
        merged_reads_mapped=Reads that were merged by bbduk and that could not map to the contigs.
    Output: 
        kmer_input=The merged fasta files.
        read_count=The total number of reads.
    """ 
    input: 
        unmerged_reads_unmapped="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_RNA/stage2/bbwrap_alignment/unmerged_reads_unmapped_joined.fasta",
        contigs="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_RNA/stage2/megahit/RNA.contigs.fa",
        merged_reads_mapped="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_RNA/stage2/bbwrap_alignment/merged_reads_unmapped.fasta"
    output: 
        kmer_input=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_RNA/stage2/kmer_input/kmer_input.fasta"),
        read_count=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_PE_RNA/stage2/kmer_input/count_kmer_input.txt")
    shell:
        """
        cat {input} > {output.kmer_input};

        echo count > {output.read_count};
        echo $(grep ">" {output.kmer_input}|wc -l) >> {output.read_count};
        """ 

rule join_unmerged_PE_DNA:
    """ 
    Rule for joining paired end reads using "N" as separator from a fastq into a fasta.
    Input: fastq=Reads that were not merged with bbmerge nor mapped to contigs from megahit.
    Output: fasta=Paired reads joined using "N" as separator.
    """ 
    input:
        fastq="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_DNA/stage1/trimming/unmerged_reads_trimmed.fq"
    output:
        fasta=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_DNA/stage2/kmer_input/unmerged_reads_joined.fasta")
    run:
        reformatted=[]
        with open(input['fastq'], 'r') as filehandle:
            contents = [line.strip() for line in filehandle.readlines()]

            for i in range(0,len(contents),8):
                headerFwd=contents[i]
                sequenceFwd=contents[i+1]

                headerRev=contents[i+4]
                sequenceRev=contents[i+5]

                readIdFwd=headerFwd.split(" ")[0]
                readIdRev=headerRev.split(" ")[0]

                readIdFwd=re.sub(r"^@", ">", readIdFwd)
                readIdRev=re.sub(r"^@", ">", readIdRev)
                
                if readIdFwd==readIdRev:
                    reformatted.append(readIdFwd)
                    reformatted.append(sequenceFwd+"N"+sequenceRev)

        with open(output['fasta'], 'w') as printresults:
            printresults.writelines("%s\n" % line for line in reformatted)            

rule convert_fasta_to_fastq_PE_DNA:
    """ 
    Rule for coverting the reads that could be merged by bbduk from fastq to fasta.
    Input: 
        fastq=file in fastq format.
    Output: 
        fasta=reads file in fasta format.
    """ 
    input: 
        fastq="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_DNA/stage1/trimming/merged_reads_trimmed.fq"
    output:
        fasta=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_DNA/stage2/kmer_input/merged_reads_trimmed.fa")
    #conda: "../../../conda/bbmap_env.yaml"
    singularity: config['singularity_bbmap_env']
    shell: 
        """
        reformat.sh \
            in={input.fastq} \
            out={output.fasta} \
            fastawrap=100000;
        """

rule create_kmer_classifier_input_PE_DNA:
    """ 
    Rule for merging reads that could not be merged by bbduk with the ones that could be merged by bbduk.
    Input:  
        unmerged_reads=Reads that could not be merged by bbduk.
        merged_reads=Reads that could be merged by bbduk.
    Output: 
        kmer_input=The merged file in fasta format.
        read_count=The number of sequences in the merged file.
    """ 
    input: 
        unmerged_reads="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_DNA/stage2/kmer_input/unmerged_reads_joined.fasta",
        merged_reads="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_DNA/stage2/kmer_input/merged_reads_trimmed.fa"
    output: 
        kmer_input=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_DNA/stage2/kmer_input/kmer_input.fasta"),
        read_count=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_PE_DNA/stage2/kmer_input/count_kmer_input.txt")
    shell:
        """
        cat {input} > {output.kmer_input};

        echo count > {output.read_count};
        echo $(grep ">" {output.kmer_input}|wc -l) >> {output.read_count};
        """ 