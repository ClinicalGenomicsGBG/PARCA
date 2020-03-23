# snakemake -rp -s main.smk --use-singularity --cores 23 --latency-wait 60
# not necessary with snakemake version 5.5.4: --singularity-args "-H /home/xerpey"
configfile: "config/config.yaml"
singularity: config['singularity_image']

shell.prefix('PATH=$PATH')


# def detect_gz_samples():
#     zipped_files=[]
#     return


rule all:
    input:
        # Quality control 
        expand("{outdir}/snakemake_results_{sample}/stage1/qc/{sample}_fastqc.html", 
                    outdir=config['outdir'],
                    sample=config['sample']),
        expand("{outdir}/snakemake_results_{sample}/stage1/qc/{sample}_fastqc.zip", 
                    outdir=config['outdir'],
                    sample=config['sample'])

include:
    "workflows/snakemake_rules/stage1_qc_trim_ec/quality_control/fastqc.smk"

# pigz -p 23 -dc a.fastq.gz > l.fq

# rule create_output_folder:
#     output:
#         "{outdir}/"
#     params:
#         output_folder_date = config['output_folder_date']
#     run:
#         if params.output_folder_date == "":
#             from datetime import date
#             today = date.today()
#             y_m_d = today.strftime("%Y_%m_%d")
#         else:
#             y_m_d = params.output_folder_date
#         shell("mkdir ")
#
# date_sample= "{date}_{sample}".format(date=y_m_d, sample=config['sample'])
# run_outdir="{outdir}/{date_sample}".format(outdir=config['outdir'],date_sample=date_sample)




