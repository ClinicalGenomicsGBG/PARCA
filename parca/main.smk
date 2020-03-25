#  snakemake -nrp -s main.smk --use-singularity --use-conda --cores 23 --latency-wait 60
# not necessary with snakemake version 5.5.4: --singularity-args "-H /home/xerpey"
# #shell.prefix('PATH=$PATH;')
import glob
from workflows.utils.setup import SetUp

configfile: "config/config.yaml"
singularity: config['singularity_image']


run = SetUp(raw_sample_dir=config['sampledir'], filenames=config['sample'])
sample = run.detect_samples()

rule all:
    input:
        expand("{outdir}/snakemake_results_{sample}/stage1/pollux/trimmed_reads.corrected.fq",
                outdir=config['outdir'], 
                sample=sample)
        ##FASTQC
        # expand("{outdir}/snakemake_results_{sample}/stage1/qc/{sample}_fastqc.{fmt}", 
        #     outdir=config['outdir'], 
        #     sample=sample, fmt=["zip","html"])
        ##TRIMMING
        # expand("{outdir}/snakemake_results_{sample}/stage1/trimming/trimmed_reads.fq", 
        #         outdir=config['outdir'], 
        #         sample=sample),
        # expand("{outdir}/snakemake_results_{sample}/stats/stage1/trimming/bbduk_stats.txt",
        #         outdir=config['outdir'], 
        #         sample=sample)

include:
    "workflows/snakemake_rules/stage1_qc_trim_ec/unzip/unzip.smk"

include:
    "workflows/snakemake_rules/stage1_qc_trim_ec/quality_control/fastqc.smk"

include:
    "workflows/snakemake_rules/stage1_qc_trim_ec/trimming/bbduk_trimming.smk"


#Negative control script


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




