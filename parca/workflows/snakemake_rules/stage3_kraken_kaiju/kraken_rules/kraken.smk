#from collections import defaultdict

rule kraken:
    input: 
        kmer_input="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage2/kmer_input/kmer_input.fasta"
    output:  
        kraken="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage3/kraken/kraken_log_{kraken_db}.txt"
    params:
        #kraken_path=config['kraken_path'],
        kraken_db_base_path=runinfo_dict['kraken_db_base_path'] #config['kraken_db_base_path']
    threads: 110
    conda: "../../../conda/kraken_kaiju_env.yaml" #config['conda_environment']
    log:"{outdir}/snakemake_results_{sample}/logs_{sample_type}_{nucleotide}/stage3/kraken/kraken_log_{kraken_db}.log"
    shell:
        """
        kraken \
            --db {params.kraken_db_base_path}/{wildcards.kraken_db} \
            --preload \
            --threads {threads} \
            --fasta-input {input.kmer_input} \
            --output {output.kraken} \
            &>{log}; 
        """
        # """
        # {params.kraken_path}/kraken \
        #     --db {params.kraken_db_base_path}/{wildcards.kraken_db} \
        #     --preload \
        #     --threads {threads} \
        #     --fasta-input {input.kmer_input} \
        #     --output {output.kraken} \
        #     &>{log}; 
        # """

rule kraken_filter_score:
    input: 
        kraken="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage3/kraken/kraken_log_{kraken_db}.txt"
    output:  
        filtered="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage3/kraken/kraken_{kraken_db}_filter_{db_limits}.txt"
    params:
        #kraken_path=config['kraken_path'],
        kraken_db_base_path=runinfo_dict['kraken_db_base_path'] #config['kraken_db_base_path']
    conda: "../../../conda/kraken_kaiju_env.yaml" #config['conda_environment']
    shell:
        """
        kraken-filter \
            --db {params.kraken_db_base_path}/{wildcards.kraken_db} \
            --threshold {wildcards.db_limits} \
            {input.kraken} \
            > {output.filtered};
        """
        # """
        # {params.kraken_path}/kraken-filter \
        #     --db {params.kraken_db_base_path}/{wildcards.kraken_db} \
        #     --threshold {wildcards.db_limits} \
        #     {input.kraken} \
        #     > {output.filtered};
        # """

rule kraken_filter_classified_RNA:
    input:
        files=expand("{{outdir}}/snakemake_results_{{sample}}/{{sample_type}}_RNA/stage3/kraken/kraken_{kraken_db}_filter_{db_limits}.txt", 
            zip,
            kraken_db=config['krakendb_RNA'], 
            db_limits=config['krakendblimits_RNA'])
    output:
        classified_filtered="{outdir}/snakemake_results_{sample}/{sample_type}_RNA/stage3/kraken/kraken_filtered_classified.txt",
        read_count="{outdir}/snakemake_results_{sample}/stats_{sample_type}_RNA/stage3/kraken/count_kraken_filtered_classified.txt"
    params:
        program="kraken"
    conda: "../../../conda/R_env.yaml" #config['conda_environment']
    script:
        "../../../scripts/kmer_processing/filter_classified.R"

rule kraken_filter_classified_DNA:
    """
    Filter best classified read.
    """
    input:
        files=expand("{{outdir}}/snakemake_results_{{sample}}/{{sample_type}}_DNA/stage3/kraken/kraken_{kraken_db}_filter_{db_limits}.txt", 
            zip,
            kraken_db=config['krakendb_DNA'], 
            db_limits=config['krakendblimits_DNA'])
    output:
        classified_filtered="{outdir}/snakemake_results_{sample}/{sample_type}_DNA/stage3/kraken/kraken_filtered_classified.txt",
        read_count="{outdir}/snakemake_results_{sample}/stats_{sample_type}_DNA/stage3/kraken/count_kraken_filtered_classified.txt"
    params:
        program="kraken"
    conda: "../../../conda/R_env.yaml" #config['conda_environment']
    script:
        "../../../scripts/kmer_processing/filter_classified.R"

# rule kraken2_SE:
#     input: 
#         kmer_input="{outdir}/snakemake_results_{sample}/SE_{nucleotide}/stage2/kmer_input/kmer_input.fasta"
#     output:  
#         kraken="{outdir}/snakemake_results_{sample}/SE_{nucleotide}/stage3/kraken2/kraken_out_{kraken_db}_{db_limits}.txt",
#         report="{outdir}/snakemake_results_{sample}/SE_{nucleotide}/stage3/kraken2/kraken_report_{kraken_db}_{db_limits}.txt"
#     params:
#         kraken2_db_base_path=config['kraken2_db_base_path']
#     threads: 110
#     conda: config['conda_environment']
#     log:"{outdir}/snakemake_results_{sample}/logs_SE_{nucleotide}/stage3/kraken2/kraken_log_{kraken_db}_{db_limits}.log"
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