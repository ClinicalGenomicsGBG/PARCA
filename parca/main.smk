import glob,re
from workflows.utils.FileProcessing import ProcessFiles
from workflows.utils.Setup import Setup

configfile: "config/config.yaml"
#snakemake -rp -s main.smk --cluster-config config/cluster.yaml --profile qsub_profile
#snakemake --dag -s main.smk| dot -Tpng > dag.png

# Read the runinfo file containg parameters for the current run.
#runinfo = ProcessFiles(config['runinfo'])
#runinfo_dict=runinfo.readYaml()

#sample_paths_dict = runinfo_dict['samplePath']
#RNA = runinfo_dict['RNA']
singularity: config['singularity_image']

# Generate settings with correct naming.
# SU=Setup(sample_paths_dict, runinfo_dict['generateSampleID'])
# settings_dict = SU.generateSettingsLists()


# print("\n\t\t~~~~~~~~ P a R C A ~~~~~~~~")
# print("\t**** Pathogen Research in Clinical Applications ****")
# print("\n**** PaRCA started for the following samples: ****")
# sample_id_list=[]
# sample_type_list=[]
# nucleotide_list=[]
# for key in settings_dict:
#     sample_id_list.append(key)
#     sample_type_list.append(settings_dict[key][1])
#     nucleotide_list.append(settings_dict[key][2])
#     print("SAMPLE ID:", key)
#     print("\tInput files:", settings_dict[key][0][0:2])
#     print("\tSample type:", settings_dict[key][1])
#     print("\tNucleotide:", settings_dict[key][2])

# print("\nResults are placed in:", runinfo_dict['outdir'], "\n")

print(expand("{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/all_classed_read_taxid_names.txt",
            zip,
            outdir=[config['outdir']]*len(sample_id_list),
            sample=sample_id_list,
            sample_type=sample_type_list,
            nucleotide=nucleotide_list
            ))
print(settings_dict)

rule all:
    input:
        # expand("{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage2/kmer_input/kmer_input.fasta",
        #     zip,
        #     outdir=[runinfo_dict['outdir']]*len(sample_id_list),
        #     sample= sample_id_list,
        #     sample_type = sample_type_list,
        #     nucleotide = nucleotide_list
        #     )
        expand("{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/all_classed_read_taxid_names.txt",
            zip,
            outdir=[config['outdir']]*len(sample_id_list),
            sample=sample_id_list,
            sample_type=sample_type_list,
            nucleotide=nucleotide_list
            )

        

##STAGE 1
include:
    "workflows/snakemake_rules/stage1_qc_trim_ec/setup/setup.smk"
include:
    "workflows/snakemake_rules/stage1_qc_trim_ec/quality_control/fastqc.smk"
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

##STAGE 3
include:
    "workflows/snakemake_rules/stage3_kraken_kaiju/kraken_rules/kraken.smk"
include:
    "workflows/snakemake_rules/stage3_kraken_kaiju/kaiju_rules/kaiju.smk"

##STAGE 4
include:
    "workflows/snakemake_rules/stage4_parse_hits/parse_hits.smk"
include:
    "workflows/snakemake_rules/stage4_parse_hits/taxonomy_processing.smk" 

##STAGE 5
include:
    "workflows/snakemake_rules/stage5_blast_processing/blast_processing.smk" 


##STAGE 6
include:
   "workflows/snakemake_rules/stage6_blast_sliced_db/blast_above_species_classed.smk"

##STAGE 7 
include:
   "workflows/snakemake_rules/stage7_blast_remaining_reads/blast_remaining.smk"

##STAGE 8
include:
    "workflows/snakemake_rules/stage8_format_results/format_results.smk"





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




