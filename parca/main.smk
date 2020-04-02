#  snakemake -nrp -s main.smk --use-singularity --use-conda --cores 23 --latency-wait 60
# not necessary with snakemake version 5.5.4: --singularity-args "-H /home/xerpey"
# snakemake --rulegraph a_b.txt | dot -Tpng > rulegraph.png
# #shell.prefix('PATH=$PATH;')
import glob
from workflows.utils.setup import SetUp

configfile: "config/config.yaml"
singularity: config['singularity_image']

def find_samples(raw_sample_dir, samples, suffix_fwd):
    if samples == "" or samples == None:
        SAMPLENAMES, = glob_wildcards(raw_sample_dir+"/{sample}"+suffix_fwd)
    else:
        SAMPLENAMES = samples
    return SAMPLENAMES

def settings(RNA="", raw_sample_dir="", samples="", suffix_fwd=".fastq", suffix_rev=""):
    try:
        RNA = RNA.capitalize() 
    except: 
        pass
    
    if RNA == "" or RNA == None or RNA == False or RNA != True or RNA=="False" or RNA=="No":
        nucleotide="DNA"
    else:
        nucleotide="RNA"
    
    if suffix_rev == "" or suffix_rev == None:
        sample_type="SE"
        sample_ids=find_samples(raw_sample_dir, samples, suffix_fwd)
    else:
        sample_type="PE"
        sample_ids=find_samples(raw_sample_dir, samples, suffix_fwd)

    return nucleotide, sample_type, sample_ids


RNA = config['RNA']
raw_sample_dir=config['sampledir']
samples=config['samples']
suffix_fwd = config['suffix_fwd']
suffix_rev = config['suffix_rev']

nucleotide, sample_type, sample_ids = settings(
            RNA,
            raw_sample_dir, 
            samples, 
            suffix_fwd, 
            suffix_rev)

rule all:
    input: 
        expand("{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage2/kmer_input/kmer_input.fasta",
            outdir=config['outdir'],
            sample=sample_ids,
            sample_type=sample_type,
            nucleotide=nucleotide),
        expand("{outdir}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage2/kmer_input/count_kmer_input.txt",
            outdir=config['outdir'],
            sample=sample_ids,
            sample_type=sample_type,
            nucleotide=nucleotide)     



# rule all:
#     input:
#         expand("{outdir}/snakemake_results_{sample}/stage2/pileup/bbmap_cov.txt",
#                 outdir=config['outdir'], 
#                 sample=sample)

        # expand("{outdir}/snakemake_results_{sample}/stage2/megahit",
        #         outdir=config['outdir'], 
        #         sample=sample),
        # expand("{outdir}/snakemake_results_{sample}/stats/stage2/megahit/assembled_contig_count.txt",
        #         outdir=config['outdir'], 
        #         sample=sample)
    
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
        ##POLLUX
        # expand("{outdir}/snakemake_results_{sample}/stage1/pollux/trimmed_reads.corrected.fq",
        #         outdir=config['outdir'], 
        #         sample=sample)
        ##FIONA
        # expand("{outdir}/snakemake_results_{sample}/stage1/fiona/trimmed_reads_fiona.fq",
        #         outdir=config['outdir'], 
        #         sample=sample)


##STAGE 1
include:
    "workflows/snakemake_rules/stage1_qc_trim_ec/setup/setup.smk"
# include:
#     "workflows/snakemake_rules/stage1_qc_trim_ec/quality_control/fastqc.smk"
include:
    "workflows/snakemake_rules/stage1_qc_trim_ec/trimming/bbduk_trimming.smk"
include:
    "workflows/snakemake_rules/stage1_qc_trim_ec/ec_pollux/ec_pollux.smk"
include:
    "workflows/snakemake_rules/stage1_qc_trim_ec/ec_fiona/ec_fiona.smk"

##STAGE 2
include:
    "workflows/snakemake_rules/stage2_assembly/megahit/megahit.smk"
include:
    "workflows/snakemake_rules/stage2_assembly/bbwrap_alignment/bbwrap_alignment.smk"
include:
    "workflows/snakemake_rules/stage2_assembly/merge_contigs_unmapped/merge_contigs_unmapped.smk"
include:
    "workflows/snakemake_rules/stage3_kraken_kaiju/kraken_rules/kraken.smk"
include:
    "workflows/snakemake_rules/stage3_kraken_kaiju/kaiju_rules/kaiju.smk"


# rule add_negative_control:
#     input: ""
#     output: ""
#     shell: ""


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




