#import glob,re

#snakemake -rp -s main.smk --cluster-config config/cluster.yaml --profile qsub_profile
#snakemake --dag -s main.smk| dot -Tpng > dag.png
# from workflows.utils.FileProcessing import ProcessFiles
# from workflows.utils.Setup import Setup

# Read the runinfo file containg parameters for the current run.
#runinfo = ProcessFiles(config['runinfo'])
#runinfo_dict=runinfo.readYaml()

#sample_paths_dict = runinfo_dict['samplePath']
#RNA = runinfo_dict['RNA']

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