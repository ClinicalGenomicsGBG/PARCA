

rule unzip_rename:
    input:
        lambda wildcards: glob.glob("{sampledir}/{sample}.f*".format(
            sampledir=config['sampledir'],
            sample=wildcards.sample
            ))
    output:
        "{outdir}/snakemake_results_{sample}/{sample}.fastq"
    threads: 23
    conda: config['conda_environment']
    shell:
        """
        if [[ {input} =~ .*\.gz$ ]]; then
            pigz -p {threads} -dc {input} > {output}
        elif [[ {input} =~ .*\.fq$ ]]; then
            ln -s {input} {output}
        elif [[ {input} =~ .*\.fastq$ ]]; then
            ln -s {input} {output}
        fi
        """