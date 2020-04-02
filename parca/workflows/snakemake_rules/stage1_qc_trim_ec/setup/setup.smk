rule unzip_rename_SE:
    input:
        lambda wildcards: "{sampledir}/{sample}{suffix_fwd}".format(
            sampledir=config['sampledir'],
            sample=wildcards.sample,
            suffix_fwd=config['suffix_fwd']
            )
    output:
        reads="{outdir}/snakemake_results_{sample}/SE_{nucleotide}/samples/{sample}.fastq",
        read_count="{outdir}/snakemake_results_{sample}/stats_SE_{nucleotide}/stage1/samples/count_raw_reads.txt"
    threads: 4
    conda: config['conda_environment']
    shell:
        """
        if [[ {input} =~ .*\.gz$ ]]; then \
            pigz -p {threads} -dc {input} > {output.reads}; \
        else \
            if [[ ! -f {output.reads} ]]; then \
                ln -s {input} {output.reads}; fi; \
        fi;  
        echo $(cat {output.reads}|wc -l)/4|bc  > {output.read_count}
        """

rule unzip_rename_PE:
    input:
        fwd=lambda wildcards: "{sampledir}/{sample}{suffix_fwd}".format(
            sampledir=config['sampledir'],
            sample=wildcards.sample,
            suffix_fwd=config['suffix_fwd']
            ),
        rev=lambda wildcards: "{sampledir}/{sample}{suffix_rev}".format(
            sampledir=config['sampledir'],
            sample=wildcards.sample,
            suffix_rev=config['suffix_rev']
            )
    output:
        fwd= "{outdir}/snakemake_results_{sample}/PE_{nucleotide}/samples/{sample}_R1.fastq",
        rev= "{outdir}/snakemake_results_{sample}/PE_{nucleotide}/samples/{sample}_R2.fastq"
    threads: 4
    conda: config['conda_environment']
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
    input: 
        fwd= "{outdir}/snakemake_results_{sample}/PE_{nucleotide}/samples/{sample}_R1.fastq",
        rev= "{outdir}/snakemake_results_{sample}/PE_{nucleotide}/samples/{sample}_R2.fastq"
    output: 
        interleaved="{outdir}/snakemake_results_{sample}/PE_{nucleotide}/samples/{sample}_interleaved.fastq",
        read_count="{outdir}/snakemake_results_{sample}/stats_PE_{nucleotide}/stage1/samples/count_raw_reads.txt"
    log: "{outdir}/snakemake_results_{sample}/logs_PE_{nucleotide}/stage1/{sample}_interleaved.log"
    benchmark: "{outdir}/snakemake_results_{sample}/benchmarks_PE_{nucleotide}/stage1/{sample}_interleaved.txt"
    conda: config['conda_environment']
    shell:
        """
        reformat.sh in={input.fwd} in2={input.rev} out={output.interleaved} > {log} 2> {log};
        echo $(cat {output.interleaved}|wc -l)/4|bc  > {output.read_count}    
        """