
import os
import yaml
import click
import snakemake
import subprocess
import pandas as pd
from workflow.utils.generate_outdir import GenerateOutdir

@click.group()
def main():
    pass
    
@main.command()
@click.option('--metadata', 'metadata',
                nargs=-1,
                required=True,
                type=click.Path(exists=True))
@click.option('--runinfo', 'runinfo',
                nargs=-1,
                required=True,
                type=click.Path(exists=True))  
@click.option('-o', '--outdir', 'outdir', type=click.Path(exists=True),
              required=True,
              help='Give a full path to a directory where all results will be placed')
@click.option('-gs', '--generate_subdir', 'generate_subdir', is_flag=True,
              help='Generate a subfolder with date within outdir a given outdir')
@click.option('-d', '--dryrun', 'dryrun', is_flag=True, help='dryrun')
def run(metadata, runinfo, dryrun, outdir, generate_subdir):
    """

    """
    if generate_subdir and outdir:
        base_outdir = outdir
        sub_outdir = GenerateOutdir.get_date_and_randomizer()
        outdir = os.path.join(base_outdir, sub_outdir)

    if not dryrun and generate_subdir and outdir:
        return_code = subprocess.call(['mkdir', outdir])
        if return_code != 0:
            raise SystemExit('Output directory could not be created')

    # config_dict_added = {
    #             'sample_paths': list(inputfiles),
    #             'threads': threads,
    #             'outdir': outdir}
    # status = snakemake.snakemake(snakefile=f'{work_dir}/Snakefile',
    #                              cores=threads,
    #                              config=config_dict_added,
    #                              workdir=work_dir,
    #                              latency_wait=30,
    #                              dryrun=dryrun)
    # print("STATUSCODE:", status)  # True or False

    #clean up if error... do not use this since if the generate outdir is not used it will remove unnecessary things 
    #return_code=subprocess.call(['rmdir', outdir])

@main.command()
def program_versions():
    """Print dependencies."""
    import subprocess
    git_version = subprocess.check_output(['git', 'describe', '--always', '--dirty'],
                                            cwd=work_dir,
                                            universal_newlines=True).strip()
    extra = ""
    if "dirty" in git_version:
        git_version = git_version.strip("-dirty")
        extra = "(working directory contains modifications)"
    print("Git repo commit id:\n\t", git_version, "\t", extra)

    # Add conda package versions

if __name__ == '__main__':
    #Create output directory from randomly generated name
    work_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Call the click groups.
    main()


#main.add_command(function)

#a=$(find $(pwd)  -type f | tr "\n" " ")
#python3 poppunk_cli.py fit-model --outdir /home/xerpey/tutorials/poppunk_tutorial/results --generate_subdir  $a --dryrun



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

# with open(self.filename, 'r') as yamlFile:
#     yamlDict = yaml.safe_load(yamlFile)

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

## dryrun functionality
## print dag functionality