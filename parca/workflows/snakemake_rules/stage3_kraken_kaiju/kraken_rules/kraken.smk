#from collections import defaultdict

rule kraken:
    """ 
    Rule for running kraken.
    Input: 
        kmer_input=A fasta file with sequences.
    Params: 
        kraken_db_base_path=Path to kraken databases.
    Output: 
        kraken=Sequence classification summary.
    """ 
    input: 
        kmer_input="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage2/kmer_input/kmer_input.fasta"
    output:  
        kraken=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage3/kraken/kraken_log_{kraken_db}.txt")
    params:
        #kraken_path=config['kraken_path'],
        kraken_path=config['kraken_path'],
        kraken_db_base_path=config['kraken_db_base_path'] #config['kraken_db_base_path']
    threads: 110
    #conda: "../../../conda/kraken_kaiju_env.yaml" #config['conda_environment']
    log:"{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/logs_{sample_type}_{nucleotide}/stage3/kraken/kraken_log_{kraken_db}.log"
    benchmark:"{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/benchmarks_{sample_type}_{nucleotide}/stage3/kraken/kraken_log_{kraken_db}.log"
    shell:
        """
        {params.kraken_path}/kraken \
            --db {params.kraken_db_base_path}/{wildcards.kraken_db} \
            --preload \
            --threads {threads} \
            --fasta-input {input.kmer_input} \
            --output {output.kraken} \
            &>{log}; 
        """
        # """
        # kraken \
        #     --db {params.kraken_db_base_path}/{wildcards.kraken_db} \
        #     --preload \
        #     --threads {threads} \
        #     --fasta-input {input.kmer_input} \
        #     --output {output.kraken} \
        #     &>{log}; 
        # """

rule kraken_filter_score:
    """ 
    Rule for filtering the kraken classifications.
    Input: Sequence classification summary.
    Params: 
        kraken_db_base_path=Path to kraken databases.
    Output:
        Filtered sequence classificaion summary. 
    """ 
    input: 
        kraken="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage3/kraken/kraken_log_{kraken_db}.txt"
    output:  
        filtered=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage3/kraken/kraken_{kraken_db}_filter_{db_limits}.txt")
    params:
        #kraken_path=config['kraken_path'],
        kraken_path=config['kraken_path'],
        kraken_db_base_path=config['kraken_db_base_path'] #config['kraken_db_base_path']
    #conda: "../../../conda/kraken_kaiju_env.yaml" #config['conda_environment']
    log: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/logs_{sample_type}_{nucleotide}/stage3/kraken/kraken_{kraken_db}_filter_{db_limits}.log"
    benchmark: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/benchmarks_{sample_type}_{nucleotide}/stage3/kraken/kraken_{kraken_db}_filter_{db_limits}.log"
    shell:
        """
        {params.kraken_path}/kraken-filter \
            --db {params.kraken_db_base_path}/{wildcards.kraken_db} \
            --threshold {wildcards.db_limits} \
            {input.kraken} \
            > {output.filtered} 2> {log};
        """
        # """
        # kraken-filter \
        #     --db {params.kraken_db_base_path}/{wildcards.kraken_db} \
        #     --threshold {wildcards.db_limits} \
        #     {input.kraken} \
        #     > {output.filtered};
        # """

rule kraken_filter_classified_RNA:
    """
    Rule for filtering all kraken classifications for the highest calculated (seq_length-kmer_len)*score+0.5.
    score = C/Q, where C is the number of k-mers mapped to LCA values in the clade rooted at the label, and Q is the number of k-mers in the sequence that lack an ambiguous nucleotide.
    Input: 
        Kraken classifications for all databases.
    Params: 
        program=Input to R-script which classifier is used.
    Output: 
        classified_filtered=Filtered Kraken classifications.
    """ 
    input:
        files=expand("{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{{sample}}/{{sample_type}}_RNA/stage3/kraken/kraken_{kraken_db}_filter_{db_limits}.txt", 
            zip,
            kraken_db=config['krakendb_RNA'], 
            db_limits=config['krakendblimits_RNA'])
    output:
        classified_filtered=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_RNA/stage3/kraken/kraken_filtered_classified.txt"),
        read_count=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_{sample_type}_RNA/stage3/kraken/count_kraken_filtered_classified.txt")
    params:
        program="kraken"
    #conda: "../../../conda/R_env.yaml" #config['conda_environment']
    log: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/logs_{sample_type}_RNA/stage3/kraken/count_kraken_filtered_classified.log"
    benchmark: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/benchmarks_{sample_type}_RNA/stage3/kraken/count_kraken_filtered_classified.log"
    singularity: config['singularity_R_env']
    script:
        "../../../scripts/kmer_processing/filter_classified.R"

rule kraken_filter_classified_DNA:
    """
    Rule for filtering all kraken classifications for the highest calculated (seq_length-kmer_len)*score+0.5.
    score = C/Q, where C is the number of k-mers mapped to LCA values in the clade rooted at the label, and Q is the number of k-mers in the sequence that lack an ambiguous nucleotide.
    Input: 
        Kaiju classifications for all databases.
    Params: 
        program=Input to R-script which classifier is used.
    Output: 
        classified_filtered=Filtered kaiju classifications.
    """ 
    input:
        files=expand("{{outdir}}/{{start_date}}_{{run_id}}/snakemake_results_{{sample}}/{{sample_type}}_DNA/stage3/kraken/kraken_{kraken_db}_filter_{db_limits}.txt", 
            zip,
            kraken_db=config['krakendb_DNA'], 
            db_limits=config['krakendblimits_DNA'])
    output:
        classified_filtered=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_DNA/stage3/kraken/kraken_filtered_classified.txt"),
        read_count=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_{sample_type}_DNA/stage3/kraken/count_kraken_filtered_classified.txt")
    params:
        program="kraken"
    #conda: "../../../conda/R_env.yaml" #config['conda_environment']
    log: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/logs_{sample_type}_DNA/stage3/kraken/count_kraken_filtered_classified.log"
    benchmark: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/benchmarks_{sample_type}_DNA/stage3/kraken/count_kraken_filtered_classified.log"
    singularity: config['singularity_R_env']
    script:
        "../../../scripts/kmer_processing/filter_classified.R"

# rule kraken2_SE:
#     input: 
#         kmer_input="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage2/kmer_input/kmer_input.fasta"
#     output:  
#         kraken="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage3/kraken2/kraken_out_{kraken_db}_{db_limits}.txt",
#         report="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage3/kraken2/kraken_report_{kraken_db}_{db_limits}.txt"
#     params:
#         kraken2_db_base_path=config['kraken2_db_base_path']
#     threads: 110
#     conda: config['conda_environment']
#     log:"{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/logs_SE_{nucleotide}/stage3/kraken2/kraken_log_{kraken_db}_{db_limits}.log"
#     shell:
#         """
#         kraken2 \
#             --confidence {wildcards.db_limits} \
#             --db {params.kraken2_db_base_path}/{wildcards.kraken_db} \
#             --threads {threads} \
#             --fasta-input {input.kmer_input} \
#             --output {output.kraken} \
#             --report {output.report} \
#             &>{log}; 
#         """ 