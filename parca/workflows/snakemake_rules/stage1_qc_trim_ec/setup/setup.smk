import pandas as pd 

rule unzip_rename_SE:
    """
    Rule for unzipping single end fastq files. If the file is unzipped it is linked to a working directory location. It is required to unzip the files since fiona error correction used in a subsequent step cannot use zipped files as input.
    Input:
        A fastq file, gzipped or not zipped. Note that other extensions e.g. ".zip" is currently not accounted for. Add if necessary.
    Output:
        An unzipped fastq file.
    """
    input:
        lambda wildcards: metadata_df.loc[metadata_df['sample_id'] == wildcards.sample].get('path_to_file').item()
        #lambda wildcards: settings_dict[wildcards.sample][0][0]
    output:
        reads=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/samples/{sample}.fastq"),
        read_count=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_SE_{nucleotide}/stage1/samples/count_raw_reads.txt"),
        read_rawpath=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_SE_{nucleotide}/stage1/samples/paths_raw_reads.txt")
    threads: 4
    #conda: config['bbmap_environment']
    log: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/logs_SE_{nucleotide}/stage1/unzip_rename_SE.log"
    # params:
    #     run_outdir='{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/samples',
    #     stats_outdir='{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_SE_{nucleotide}/stage1/samples'
    shell:
        """
        if [[ {input} =~ .*\.fastq\.gz$ || {input} =~ .*\.fq\.gz$ ]]; then \
            pigz -p {threads} -dc {input} > {output.reads} 2> {log}; \
        elif [[ {input} =~ .*\.fastq$ || {input} =~ .*\.fq$ ]]; then \
            ln -s {input} {output.reads} 2> {log}; \
        else \
            echo "Unsupported file extension: Input should have extension .fastq or .fq or .gz" > {log}; \
        fi;  
        echo count > {output.read_count};
        echo $(cat {output.reads}|wc -l)/4|bc  >> {output.read_count} 2> {log};
        
        echo type > {output.read_rawpath};
        echo {input} >> {output.read_rawpath};
        """
        # """
        # if [[ {input} =~ .*\.gz$ ]]; then \
        #     pigz -p {threads} -dc {input} > {output.reads}; \
        # else \
        #     if [[ ! -f {output.reads} ]]; then \
        #         ln -s {input} {output.reads}; fi; \
        # fi;  
        # echo count > {output.read_count};
        # echo $(cat {output.reads}|wc -l)/4|bc  >> {output.read_count};
        # """

rule unzip_rename_PE:
    """
    Rule for unzipping paired end end fastq files and renaming them to the same convention. If the file is unzipped it is linked to a working directory location. It is required to unzip the files since fiona error correction used in a subsequent step cannot use zipped files as input.
    Input:
        Two fastq files, gzipped or not zipped. Note that other extensions e.g. ".zip" is currently not accounted for. Add if necessary.
    Output:
        Forward and reverse unzipped fastq files.
    """
    input:
        fwd = lambda wildcards: metadata_df.loc[(metadata_df['sample_id'] == wildcards.sample) & (metadata_df['fwd_or_rev'] == 'fwd')].get('path_to_file').item(),
        rev = lambda wildcards: metadata_df.loc[(metadata_df['sample_id'] == wildcards.sample) & (metadata_df['fwd_or_rev'] == 'rev')].get('path_to_file').item()
        # lambda wildcards: settings_dict[wildcards.sample][0][0],
        # lambda wildcards: settings_dict[wildcards.sample][0][1]
    output:
        fwd = temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/samples/{sample}_R1.fastq"),
        rev = temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/samples/{sample}_R2.fastq"),
        ,
        read_rawpath=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_PE_{nucleotide}/stage1/samples/paths_raw_reads.txt")
        #read_count="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_SE_{nucleotide}/stage1/samples/count_raw_reads.txt"
    threads: 4
    #conda: config['bbmap_environment']
    log: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/logs_PE_{nucleotide}/stage1/unzip_rename_PE.log"
    shell:
        """
        if [[ {input.fwd} =~ .*\.fastq\.gz$ || {input.fwd} =~ .*\.fq\.gz$ ]]; then \
            pigz -p {threads} -dc {input.fwd} > {output.fwd} 2> {log}; \
        elif [[ {input.fwd} =~ .*\.fastq$ || {input.fwd} =~ .*\.fq$ ]]; then \
            ln -s {input.fwd} {output.fwd} 2> {log}; \
        else \
            echo "Unsupported file extension: Input should have extension .fastq or .fq or .gz"; \
        fi;  


        if [[ {input.rev} =~ .*\.fastq\.gz$ || {input.rev} =~ .*\.fq\.gz$ ]]; then \
            pigz -p {threads} -dc {input.rev} > {output.rev} 2>> {log}; \
        elif [[ {input.rev} =~ .*\.fastq$ || {input.rev} =~ .*\.fq$ ]]; then \
            ln -s {input.rev} {output.rev} 2>> {log}; \
        else \
            echo "Unsupported file extension: Input should have extension .fastq or .fq or .gz"; \
        fi; 

        echo type > {output.read_rawpath};
        echo {input.fwd} >> {output.read_rawpath};
        echo {input.rev} >> {output.read_rawpath};
        """
        # """
        # if [[ {input.fwd} =~ .*\.gz$ ]]; then \
        #     pigz -p {threads} -dc {input.fwd} > {output.fwd}; \
        # else \
        #     if [[ ! -f {output.fwd} ]]; then \
        #         ln -s {input.fwd} {output.fwd}; fi; \
        # fi;  


        # if [[ {input.rev} =~ .*\.gz$ ]]; then \
        #     pigz -p {threads} -dc {input.rev} > {output.rev}; \
        # else \
        #     if [[ ! -f {output.rev} ]]; then \
        #         ln -s {input.rev} {output.rev}; fi; \
        # fi;  
        # """

rule interleave_PE:
    """
    Rule for interleving paired end files for compatibility with following rules, might be able to remove this part in following updates.
    Input:
        Unzipped forward and revese fastq-files
    Output:
        An unzipped fastq file.
    """
    input: 
        fwd= "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/samples/{sample}_R1.fastq",
        rev= "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/samples/{sample}_R2.fastq"
    output: 
        interleaved=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/samples/{sample}_interleaved.fastq"),
        read_count=temp("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_PE_{nucleotide}/stage1/samples/count_raw_reads.txt")
    log: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/logs_PE_{nucleotide}/stage1/{sample}_interleaved.log"
    benchmark: "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/benchmarks_PE_{nucleotide}/stage1/{sample}_interleaved.txt"
    conda: "../../../conda/bbmap_env.yaml"
    shell:
        """
        reformat.sh in={input.fwd} in2={input.rev} out={output.interleaved} &> {log}
        echo count > {output.read_count};
        echo $(cat {output.interleaved}|wc -l)/4|bc  >> {output.read_count};    
        """