

rule compare_kmer_results:
    input:
        kraken="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage3/kraken/kraken_filtered_classified.txt",
        kaiju="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage3/kaiju/kaiju_filtered_classified.txt"
    output:
        kraken_doublets="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/kraken_doublets.txt",
        kaiju_doublets="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/kaiju_doublets.txt",
        singletons="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/singletons.txt",
        merged="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/merged_total.txt",
        read_count_doublet="{outdir}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage4/count_doublets.txt",
        read_count_singletons="{outdir}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage4/count_singletons.txt"
    conda: config['conda_environment']
    script:
        "../../scripts/kmer_processing/compare_outputs.R"

rule merge_doublets:
    input:
        kaiju_doublets="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/kaiju_doublets.txt",
        kraken_doublets="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/kraken_doublets.txt"
    output:
        combined="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/combined_kraken_kaiju.txt"
    params:
        nodes=config['nodes']
    conda: config['conda_environment']
    log: "{outdir}/snakemake_results_{sample}/logs_{sample_type}_{nucleotide}/stage4/mergeOutputs.log"
    shell:
        """
        kaiju-mergeOutputs \
        -i {input.kaiju_doublets} \
        -j {input.kraken_doublets} \
        -o {output.combined} \
        -v -c lowest \
        -t {params.nodes} &> {log};
        """
#"mergeOutputs -i $outdir/kaiju_compare.txt -j $outdir/kraken_compare.txt -o $outdir/combined_kraken_kaiju.txt -v -c lowest -t /tmp/pathfinder_dbs/nodes.dmp";

rule add_taxon_names_doublets:
    input:
        combined="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/combined_kraken_kaiju.txt"
    output:
        named="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/combined_kraken_kaiju_names_unfiltered.txt"
    params:
        nodes=config['nodes'],
        names=config['names']
    conda: config['conda_environment']
    shell:
        """
        kaiju-addTaxonNames \
        -t {params.nodes} \
        -n {params.names} \
        -i {input.combined} \
        -o {output.named} \
        -r superkingdom,class,order,family,genus,species;
        """
