rule kraken_RNA:
    input: 
        ""
    output:
        kmer_input="{outdir}/snakemake_results_{sample}/{sample_type}_RNA/stage2/kmer_input/kmer_input.fasta",
        databases=""
    shell:
        """
        """

# @krakendb=("eukaryotes_cds","viruses","bacteria_progenomes_cds_total");
# @krakendbnames=("eukaryotes_cds","viruses","bacteria_progenomes_cds_total");
# @krakendblimits=("0.15","0.05","0.05");	

rule kraken_DNA:
    input: 
        ""
    output:
        kmer_input="{outdir}/snakemake_results_{sample}/{sample_type}_DNA/stage2/kmer_input/kmer_input.fasta",
        databases=""
    shell:
        """
        """

# @krakendb=("eukaryotes_768","viruses","bacteria_progenomes_total_768");
# @krakendbnames=("eukaryotes_768","viruses","bacteria_progenomes_total_768");
# @krakendblimits=("0.15","0.05","0.05");
