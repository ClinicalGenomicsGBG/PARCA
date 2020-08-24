
rule unzip_rename_SE:
    """
    Rule for unzipping single end fastq files. If the file is unzipped it is linked to a working directory location.
    Input:
        A fastq file, gzipped or not zipped. Note that other extensions e.g. ".zip" is currently not accounted for. Add if necessary.
    Output:
        An unzipped fastq file.
    """
    input:
        lambda wildcards: settings_dict[wildcards.sample][0]
    output:
        reads="{outdir}/snakemake_results_{sample}/SE_{nucleotide}/samples/{sample}.fastq",
        read_count="{outdir}/snakemake_results_{sample}/stats_SE_{nucleotide}/stage1/samples/count_raw_reads.txt"
    threads: 4
    #conda: config['bbmap_environment']
    shell:
        """
        if [[ {input} =~ .*\.gz$ ]]; then \
            pigz -p {threads} -dc {input} > {output.reads}; \
        else \
            if [[ ! -f {output.reads} ]]; then \
                ln -s {input} {output.reads}; fi; \
        fi;  
        echo count > {output.read_count};
        echo $(cat {output.reads}|wc -l)/4|bc  >> {output.read_count};
        """

rule unzip_rename_PE:
    """
    Rule for unzipping paired end end fastq files and renaming them to the same convention. If the file is unzipped it is linked to a working directory location.
    Input:
        Two fastq files, gzipped or not zipped. Note that other extensions e.g. ".zip" is currently not accounted for. Add if necessary.
    Output:
        Forward and reverse unzipped fastq files.
    """
    input:
        fwd=lambda wildcards: settings_dict[wildcards.sample][0][0],
        rev=lambda wildcards: settings_dict[wildcards.sample][0][1]
    output:
        fwd= "{outdir}/snakemake_results_{sample}/PE_{nucleotide}/samples/{sample}_R1.fastq",
        rev= "{outdir}/snakemake_results_{sample}/PE_{nucleotide}/samples/{sample}_R2.fastq"
    threads: 4
    #conda: config['bbmap_environment']
    shell:
        """
        if [[ {input.fwd} =~ .*\.gz$ ]]; then \
            pigz -p {threads} -dc {input.fwd} > {output.fwd}; \
        else \
            if [[ ! -f {output.fwd} ]]; then \
                ln -s {input.fwd} {output.fwd}; fi; \
        fi;  
        if [[ {input.rev} =~ .*\.gz$ ]]; then \
            pigz -p {threads} -dc {input.rev} > {output.rev}; \
        else \
            if [[ ! -f {output.rev} ]]; then \
                ln -s {input.rev} {output.rev}; fi; \
        fi;  
        """

rule interleave_PE:
    """
    Rule for interleving paired end files for compatibility with following rules, might be able to remove this part in following updates.
    Input:
        Unzipped forward and revese fastq-files
    Output:
        An unzipped fastq file.
    """
    input: 
        fwd= "{outdir}/snakemake_results_{sample}/PE_{nucleotide}/samples/{sample}_R1.fastq",
        rev= "{outdir}/snakemake_results_{sample}/PE_{nucleotide}/samples/{sample}_R2.fastq"
    output: 
        interleaved="{outdir}/snakemake_results_{sample}/PE_{nucleotide}/samples/{sample}_interleaved.fastq",
        read_count="{outdir}/snakemake_results_{sample}/stats_PE_{nucleotide}/stage1/samples/count_raw_reads.txt"
    log: "{outdir}/snakemake_results_{sample}/logs_PE_{nucleotide}/stage1/{sample}_interleaved.log"
    benchmark: "{outdir}/snakemake_results_{sample}/benchmarks_PE_{nucleotide}/stage1/{sample}_interleaved.txt"
    conda: "../../../conda/bbmap_env.yaml"
    shell:
        """
        reformat.sh in={input.fwd} in2={input.rev} out={output.interleaved} &> {log}
        echo count > {output.read_count};
        echo $(cat {output.interleaved}|wc -l)/4|bc  >> {output.read_count};    
        """