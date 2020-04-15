rule kaiju:
    input: 
        kmer_input="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage2/kmer_input/kmer_input.fasta"
    output:
        kaiju="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage3/kaiju/kaijuresults_{kaiju_db}_{kaiju_score}_{kaiju_matches}.txt"
    params:
        kaiju_db_base_path=config['kaiju_db_base_path'],
        kaijunames=config['kaiju_names']
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

rule kaiju_filter_classified_RNA:
    input:
        files=expand("{{outdir}}/snakemake_results_{{sample}}/{{sample_type}}_RNA/stage3/kaiju/kaijuresults_{kaiju_db}_{kaiju_score}_{kaiju_matches}.txt",
            zip,
            kaiju_db=config['kaijudb_RNA'], 
            kaiju_score=config['kaijuscore_RNA'],
            kaiju_matches=config['kaijumatches_RNA'] ) 
    output:
        classified_filtered="{outdir}/snakemake_results_{sample}/{sample_type}_RNA/stage3/kaiju/kaiju_filtered_classified.txt",
        read_count="{outdir}/snakemake_results_{sample}/stats_{sample_type}_RNA/stage3/kaiju/count_kaiju_filtered_classified.txt"
    params:
        program="kaiju"
    conda: config['conda_environment']
    script:
            "../../../scripts/kmer_processing/filter_classified.R"

# Kaiju 1.5.0
# Copyright 2015-2017 Peter Menzel, Anders Krogh
# License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>

# Usage:
#    kaiju -t nodes.dmp -f kaiju_db.fmi -i reads.fastq [-j reads2.fastq]

# Mandatory arguments:
#    -t FILENAME   Name of nodes.dmp file
#    -f FILENAME   Name of database (.fmi) file
#    -i FILENAME   Name of input file containing reads in FASTA or FASTQ format

# Optional arguments:
#    -j FILENAME   Name of second input file for paired-end reads
#**    -o FILENAME   Name of output file. If not specified, output will be printed to STDOUT
#**    -z INT        Number of parallel threads (default: 1)
#**    -a STRING     Run mode, either "mem"  or "greedy" (default: mem)
#**    -e INT        Number of mismatches allowed in Greedy mode (default: 0)
#**    -m INT        Minimum match length (default: 11)
#**    -s INT        Minimum match score in Greedy mode (default: 65)
#    -E FLOAT      Minimum E-value in Greedy mode
#**    -x            Enable SEG low complexity filter
#    -p            Input sequences are protein sequences
#**    -v            Enable verbose output            
