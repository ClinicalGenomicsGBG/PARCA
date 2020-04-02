rule create_kmer_classifier_input_SE_RNA:
    input:
        contigs="{outdir}/snakemake_results_{sample}/SE_RNA/stage2/megahit/RNA.contigs.fa",
        unmapped_reads="{outdir}/snakemake_results_{sample}/SE_RNA/stage2/bbwrap_alignment/unmapped_reads.fasta"
    output:
        kmer_input="{outdir}/snakemake_results_{sample}/SE_RNA/stage2/kmer_input/kmer_input.fasta",
        read_count="{outdir}/snakemake_results_{sample}/stats_SE_RNA/stage2/kmer_input/count_kmer_input.txt"
    conda: config['conda_environment']
    shell:
        """
        cat {input.contigs} {input.unmapped_reads} > {output.kmer_input};
        echo $(cat {output.kmer_input}|wc -l)/4|bc  > {output.read_count};
        """
