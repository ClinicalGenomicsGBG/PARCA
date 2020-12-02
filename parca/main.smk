#import glob,re
import yaml
import pandas as pd
import numpy as np

configfile: "config/config.yaml"

singularity: config['singularity_image']

run_dict = config['run_dict']
metadata_dict = config['metadata_dict']

metadata_dataframe = pd.DataFrame(metadata_dict)

outdir = config['outdir']

def generate_pipeline_input(list_run_dictionary):
    pipeline_input=[]
    for run in list_run_dictionary:
        if "case" in run and "control" in run:
            case_sample_id = run['case']
            control_sample_id = run['control']
            case_control_list=expand(f'{outdir}/case_{case_sample_id}_control_{control_sample_id}')
            pipeline_input.append(case_control_list)
        elif "case" in run:
            case_sample_id = run['case']
            case_list=expand(f'{outdir}/case_{case_sample_id}')
            pipeline_input.append(case_list)
    return pipeline_input

rule all:
    run:
        print(generate_pipeline_input(run_dict))


rule control_and_case:
    input:
        case=expand("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/all_classed_read_taxid_names.txt", 
            outdir="", sample="", sample_type="", nucleotide="")
    output:
        "{outdir}/{start_date}_{run_id}/krona.txt"
    shell:
        """
        touch {output}
        """

rule case:
    output:
        "{outdir}/case_{case}.txt"
    shell:
        """
        touch {output}
        """

# Rule all
print(expand("{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/all_classed_read_taxid_names.txt",
#             zip,
#             outdir=[config['outdir']]*len(sample_id_list),
#             sample=sample_id_list,
#             sample_type=sample_type_list,
#             nucleotide=nucleotide_list
#             ))


# rule all:
#     input:
#         # expand("{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage2/kmer_input/kmer_input.fasta",
#         #     zip,
#         #     outdir=[config['outdir']]*len(sample_id_list),
#         #     sample= sample_id_list,
#         #     sample_type = sample_type_list,
#         #     nucleotide = nucleotide_list
#         #     )
#         expand("{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/all_classed_read_taxid_names.txt",
#             zip,
#             outdir=[config['outdir']]*len(sample_id_list),
#             sample=sample_id_list,
#             sample_type=sample_type_list,
#             nucleotide=nucleotide_list
#             )

        

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


#run_dict = [{'run_id': 'run_1', 'case': 'sample_1', 'control': 'sample_2'}, {'run_id': 'run_2', 'case': 'sample_2'}]
#metadata_dict = [{'sample_id': 'sample_1', 'start_date': 20201104, 'nucleotide': 'RNA', 'fwd_or_rev': 'fwd', 'path_to_file': '/apps/bio/dev_repos/parca/demo/raw_samples/SRR1761912_1.fastq.gz', 'adapters': np.nan, 'PE_or_SE': 'PE'}, {'sample_id': 'sample_1', 'start_date': 20201104, 'nucleotide': 'RNA', 'fwd_or_rev': 'rev', 'path_to_file': '/apps/bio/dev_repos/parca/demo/raw_samples/SRR1761912_2.fastq.gz', 'adapters': np.nan, 'PE_or_SE': 'PE'}, {'sample_id': 'sample_2', 'start_date': 20201104, 'nucleotide': 'DNA', 'fwd_or_rev': 'fwd', 'path_to_file': '/apps/bio/dev_repos/parca/demo/raw_samples/a.fastq.gz', 'adapters': np.nan, 'PE_or_SE': 'SE'}]