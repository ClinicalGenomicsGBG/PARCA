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
