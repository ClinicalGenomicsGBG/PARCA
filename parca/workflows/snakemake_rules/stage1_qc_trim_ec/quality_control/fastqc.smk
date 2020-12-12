

rule quality_control_raw_SE:
    input:
        fastq="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/samples/{sample}.fastq"
    output:
        html="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage1/qc_raw/{sample}_fastqc.html",
        out_zip="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage1/qc_raw/{sample}_fastqc.zip"
    params:
        dir="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage1/qc_raw/"
    threads: 4
    conda: "../../../conda/bbmap_env.yaml"
    shell:
        """
        fastqc -o {params.dir} \
               -t {threads} \
               {input.fastq};
        """

rule quality_control_raw_PE:
    input:
        fwd="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/samples/{sample}_R1.fastq",
        rev="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/samples/{sample}_R2.fastq"
    output:
        fwd_html="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/qc_raw/{sample}_R1_fastqc.html",
        fwd_zip="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/qc_raw/{sample}_R1_fastqc.zip",
        rev_html="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/qc_raw/{sample}_R2_fastqc.html",
        rev_zip="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/qc_raw/{sample}_R2_fastqc.zip"
    params:
        dir="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/qc_raw/"
    threads: 4
    conda: "../../../conda/bbmap_env.yaml"
    shell:
        """
        fastqc -o {params.dir} \
               -t {threads} \
               {input};
        """