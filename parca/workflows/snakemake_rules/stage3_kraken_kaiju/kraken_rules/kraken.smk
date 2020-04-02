from collections import defaultdict

rule kraken:
    input: 
        kmer_input="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage2/kmer_input/kmer_input.fasta"
    output:  
        kraken="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage3/kraken/kraken_log_{kraken_db}.txt"
        #database="{outdir}/snakemake_results_{sample}/stats_{sample_type}_RNA/stage3/{}_settings.txt"
    params:
        kraken_path=config['kraken_path'],
        kraken_db_base_path=config['kraken_db_base_path']
    threads: 110
    conda: config['conda_environment']
    shell:
        """
        {params.kraken_path}/kraken \
            --db {params.kraken_db_base_path}/{wildcards.kraken_db} \
            --preload \
            --threads {threads} \
            --fasta-input {input.kmer_input} \
            --output {output.kraken}"; 
        """

rule kraken_filter:
    input: 
        kraken="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage3/kraken/kraken_log_{kraken_db}.txt"
    output:  
        filtered="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage3/kraken/kraken_{kraken_db}_filter_{db_limits}.txt"
    params:
        kraken_path=config['kraken_path'],
        kraken_db_base_path=config['kraken_db_base_path']
    threads: 110
    conda: config['conda_environment']
    shell:
        """
        {params.kraken_path}/kraken-filter \
            --db {params.kraken_db_base_path}/{wildcards.kraken_db} \
            --threshold {wildcards.db_limits} \
            {input.kraken} \
            > {output.filtered}";
        """

rule filter_kraken_classified_kraken_RNA:
    input:
        krakenfile=expand("{{outdir}}/snakemake_results_{{sample}}/{{sample_type}}_RNA/stage3/kraken/kraken_{kraken_db}_filter_{db_limits}.txt", 
            zip,
            kraken_db=config['krakendb_RNA'], 
            db_limits=config['krakendblimits_RNA'])
    output:
        kraken_classified_filtered="{outdir}/snakemake_results_{sample}/{sample_type}_RNA/stage3/kraken/kraken_filtered_classified.txt",
        read_count="{outdir}/snakemake_results_{sample}/stats_{sample_type}_RNA/stage3/kraken/count_kraken_filtered_classified.txt"
    params:
        kmer_len=30
    conda: config['conda_environment']
    script:
        """
        ../../../scripts/kmer_processing/kraken_filter_classified.R
        """

rule filter_kraken_classified_kraken_DNA:
    input:
        krakenfile=expand("{{outdir}}/snakemake_results_{{sample}}/{{sample_type}}_DNA/stage3/kraken/kraken_{kraken_db}_filter_{db_limits}.txt", 
            zip,
            kraken_db=config['krakendb_DNA'], 
            db_limits=config['krakendblimits_DNA'])
    output:
        kraken_classified_filtered="{outdir}/snakemake_results_{sample}/{sample_type}_DNA/stage3/kraken/kraken_filtered_classified.txt",
        read_count="{outdir}/snakemake_results_{sample}/stats_{sample_type}_DNA/stage3/kraken/count_kraken_filtered_classified.txt"
    params:
        kmer_len=30
    conda: config['conda_environment']
    script:
        """
        ../../../scripts/kmer_processing/kraken_filter_classified.R
        """



# foreach $currentdb (@krakendb){
#     ($acceptedafterfilter, $scorehashref)=krakenprocess($currentdb, $outdir, \%krakenhash);
#     %krakenhash=%$scorehashref;
#     print "Accepted $acceptedafterfilter reads after P-value filtering from kraken $currentdb\n";
#     print STATSFILE "KrakenAfterFilter: $currentdb, $acceptedafterfilter\n";
#     $statisticshash{$sampletype}{"Kraken Filtered Reads"}{$currentdb}=$acceptedafterfilter;
#   }
# ($finalkrakenhashref)=krakencompare(\%krakenhash);
# %finalkrakenhash=%$finalkrakenhashref;
# @max=keys %finalkrakenhash;
# print "Reads in finalkrakenhash ".scalar @max."\n";
# $totalaccepted=scalar @max;
# print STATSFILE "ReadsClassedByKraken: $totalaccepted\n";
# $statisticshash{$sampletype}{"Kraken Filtered Reads"}{"Total"}=$totalaccepted;
