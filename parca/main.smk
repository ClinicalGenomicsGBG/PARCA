
import glob
import re
from workflows.utils.FileProcessing import ProcessFiles

configfile: "config/config.yaml"
singularity: config['singularity_image']


runinfo_dict = ProcessFiles(config['runinfo'])
sample_paths_dict = runinfo_dict['samplePath']

RNA = runinfo_dict['RNA']

SU=Setup(sample_paths_dict, RNA, False)
sample_id_list, sample_type_list, nucleotide_list= SU.generateSettingsLists()


rule all:
    input:
        expand("{outdir}/{}/{}",
            zip,
            outdir="",
            sample_id_list,
            sample_type_list,
            nucleotide_list
            )
        # expand("{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/all_classed_read_taxid_names.txt",
        #     outdir=config['outdir'],
        #     sample=sample_ids,
        #     sample_type=sample_type,
        #     nucleotide=nucleotide
        #     )
        
        # expand("{sampledir}/{sample}{suffix_fwd}",
        #     sampledir=raw_sample_dir,
        #     sample=sample_ids,
        #     suffix_fwd=suffix_fwd
        #     )




        #expand("{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/combined_doublets_singletons.txt",
            # outdir=config['outdir'],
            # sample=sample_ids,
            # sample_type=sample_type,
            # nucleotide=nucleotide
            # )

# ##STAGE 1
# include:
#     "workflows/snakemake_rules/stage1_qc_trim_ec/setup/setup.smk"
# # include:
# #     "workflows/snakemake_rules/stage1_qc_trim_ec/quality_control/fastqc.smk"
# include:
#     "workflows/snakemake_rules/stage1_qc_trim_ec/trimming/bbduk_trimming.smk"
# include:
#     "workflows/snakemake_rules/stage1_qc_trim_ec/ec_pollux/ec_pollux.smk"
# include:
#     "workflows/snakemake_rules/stage1_qc_trim_ec/ec_fiona/ec_fiona.smk"

# ##STAGE 2
# include:
#     "workflows/snakemake_rules/stage2_assembly/megahit/megahit.smk"
# include:
#     "workflows/snakemake_rules/stage2_assembly/bbwrap_alignment/bbwrap_alignment.smk"
# include:
#     "workflows/snakemake_rules/stage2_assembly/merge_contigs_unmapped/merge_contigs_unmapped.smk"

# ##STAGE 3
# include:
#     "workflows/snakemake_rules/stage3_kraken_kaiju/kraken_rules/kraken.smk"
# include:
#     "workflows/snakemake_rules/stage3_kraken_kaiju/kaiju_rules/kaiju.smk"

# ##STAGE 4
# include:
#     "workflows/snakemake_rules/stage4_parse_hits/parse_hits.smk"
# include:
#     "workflows/snakemake_rules/stage4_parse_hits/taxonomy_processing.smk" 

# ##STAGE 5
# include:
#     "workflows/snakemake_rules/stage5_blast_processing/blast_processing.smk" 


# ##STAGE 6
# include:
#    "workflows/snakemake_rules/stage6_blast_sliced_db/blast_above_species_classed.smk"

# ##STAGE 7 
# include:
#    "workflows/snakemake_rules/stage7_blast_remaining_reads/blast_remaining.smk"

# ##STAGE 8
# include:
#     "workflows/snakemake_rules/stage8_format_results/format_results.smk"





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




