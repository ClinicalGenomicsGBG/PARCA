

rule quality_control_raw:
    input:
        lambda wildcards: "{sampledir}/{sample}.fastq".format(
                            sampledir=config['sampledir'], 
                            sample=wildcards.sample)
    output:
        "{outdir}/snakemake_results_{sample}/stage1/qc/{sample}_fastqc.html",
        "{outdir}/snakemake_results_{sample}/stage1/qc/{sample}_fastqc.zip"
    params:
        dir="{outdir}/snakemake_results_{sample}/stage1/qc"
    # singularity: 
    #     config['singularity_image']
    threads: 23
    shell:
        """
        fastqc -o {params.dir} \
               -t {threads} \
               {input}
        """