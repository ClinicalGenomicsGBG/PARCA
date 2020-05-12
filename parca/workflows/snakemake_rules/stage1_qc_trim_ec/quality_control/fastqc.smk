

rule quality_control_raw:
    input:
        "{outdir}/snakemake_results_{sample}/{sample}.fastq"
    output:
        "{outdir}/snakemake_results_{sample}/stage1/qc/{sample}_fastqc.html",
        "{outdir}/snakemake_results_{sample}/stage1/qc/{sample}_fastqc.zip"
    params:
        dir="{outdir}/snakemake_results_{sample}/stage1/qc"
    threads: 23
    conda: config['conda_environment']
    shell:
        """
        fastqc -o {params.dir} \
               -t {threads} \
               {input};
        """