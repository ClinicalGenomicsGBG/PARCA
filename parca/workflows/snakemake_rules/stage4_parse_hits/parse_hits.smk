

rule compare_kmer_results:
    """ 
    Rule for comparing the kaiju and kraken results. 
    Input: 
        kraken=Kraken classifications.
        kaiju=Kaiju classifications.
    Output: 
        kraken_doublets=Kraken classifications for sequences that Kaiju also was able to classify.
        Kaiju_doublets=Kaiju classifications for sequences that Kraken also was able to classify.
        singletons=Classifications that only Kraken or Kaiju was able to make.
        merged=Kraken and Kaiju classifications in the same dataframe.
        read_count_doublet=The number of sequences that both Kaiju and Kraken was able to classify.
        read_count_singletons=The number of sequences that only Kaiju or Kraken was able to classify.
    """ 
    input:
        kraken="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage3/kraken/kraken_filtered_classified.txt",
        kaiju="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage3/kaiju/kaiju_filtered_classified.txt"
    output:
        kraken_doublets=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/kraken_doublets.txt"),
        kaiju_doublets=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/kaiju_doublets.txt"),
        singletons=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/singletons.txt"),
        merged=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/merged_total.txt"),
        read_count_doublet=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage4/count_doublets.txt"),
        read_count_singletons=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage4/count_singletons.txt")
    #conda: "../../conda/R_env.yaml" #config['conda_environment']
    log: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/logs_{sample_type}_{nucleotide}/stage4/compare_kmer_results.log"
    benchmark: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/benchmarks_{sample_type}_{nucleotide}/stage4/compare_kmer_results.log"
    singularity: config['singularity_R_env']
    script:
        "../../scripts/kmer_processing/compare_outputs.R"

rule merge_doublets:
    """ 
    Rule for comparing the classifications from Kraken and Kaiju so that conflicting classifications for reads are solved by taking the lowest taxonomic ID if the classifications belong to the same lineage, if they do not belong to the same lineage the LCA is selected.
    Input: 
        kaiju_doublets=Kaiju classifications for sequences that Kraken also was able to classify.
        kraken_doublets=Kraken classifications for sequences that Kaiju also was able to classify.
    Params: 
        names_nodes_dmp_dir=Names and nodes dmp file.
    Output: 
        Merged classifications from Kraken and Kaiju.
    """ 
    input:
        kaiju_doublets="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/kaiju_doublets.txt",
        kraken_doublets="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/kraken_doublets.txt"
    output:
        combined=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/combined_kraken_kaiju.txt")
    params:
        names_nodes_dmp_dir=config['names_nodes_dmp_dir'] #config['names_nodes_dmp_dir']
    #conda: "../../conda/kaiju_env.yaml" #config['conda_environment']
    singularity: config['singularity_kaiju_env']
    log: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/logs_{sample_type}_{nucleotide}/stage4/mergeOutputs.log"
    benchmark: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/benchmarks_{sample_type}_{nucleotide}/stage4/mergeOutputs.log"
    shell:
        """
        kaiju-mergeOutputs \
            -i {input.kaiju_doublets} \
            -j {input.kraken_doublets} \
            -o {output.combined} \
            -v -c lowest \
            -t {params.names_nodes_dmp_dir}/nodes.dmp &> {log};
        """
#"mergeOutputs -i $outdir/kaiju_compare.txt -j $outdir/kraken_compare.txt -o $outdir/combined_kraken_kaiju.txt -v -c lowest -t /tmp/pathfinder_dbs/nodes.dmp";


