rule kaiju:
    input: 
        kmer_input="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage2/kmer_input/kmer_input.fasta"
    output:
        kaiju="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage3/kaiju/kaijuresults_{kaiju_db}_{kaiju_score}_{kaiju_matches}.txt"
    params:
        #kaiju_path=config['kaiju_path'],
        kaiju_db_base_path=config['kaiju_db_base_path'],
        kaijunames=config['kaijunames']
    threads: 110
    conda: config['conda_environment']
    shell:
        """
        kaiju \
            -t {params.kaijunames} \
            -f {params.kaiju_db_base_path}/{wildcards.kaiju_db}.fmi \
            -i {input.kmer_input} \
            -x -v \
            -z {threads} \
            -s {wildcards.kaiju_score} \
            -m {wildcards.kaiju_matches} \
            -e 5 \
            -a greedy \
            -o {output.kaiju};
        """
