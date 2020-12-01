
import os
import click
import snakemake
import subprocess
import pandas as pd
from workflows.utils.generate_outdir import GenerateOutdir
from workflows.utils.process_runinfo_metadata import ProcessRuninfoMetadata

@click.group()
def main():
    """
    ~~~~~~~~ P a R C A ~~~~~~~~
    **** Pathogen Research in Clinical Applications ****
    **** PaRCA started for the following samples: ****
    """
    pass

@main.command()
@click.option('-m', '--metadata', 'metadata',
              nargs=-1,
              required=True,
              type=click.Path(exists=True))
@click.option('-r', '--runinfo', 'runinfo',
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
    """a"""
    if generate_subdir and outdir:
        base_outdir = outdir
        sub_outdir = GenerateOutdir.get_date_and_randomizer()
        outdir = os.path.join(base_outdir, sub_outdir)

    if not dryrun and generate_subdir and outdir:
        return_code = subprocess.call(['mkdir', outdir])
        if return_code != 0:
            raise SystemExit('Output directory could not be created')

    run_dict = ProcessRuninfoMetadata.generate_runinfo_dict(runinfo)
    metadata_dict = ProcessRuninfoMetadata.generate_metadata_dict(metadata)

    config_dict_added = {
                'run_dict': run_dict,
                'metadata_dict': metadata_dict,
                'outdir': outdir}

    status = snakemake.snakemake(snakefile=f'{work_dir}/Snakefile',
                                 config=config_dict_added,
                                 workdir=work_dir,
                                 latency_wait=30,
                                 dryrun=dryrun)

    # print("STATUSCODE:", status)  # True or False

    #clean up if error... do not use this since if the generate outdir is not used it will remove unnecessary things 
    #return_code=subprocess.call(['rmdir', outdir])


@main.command()
def program_versions():
    """Print dependencies."""
    import subprocess
    git_version = subprocess.check_output(['git', 'describe',
                                          '--always', '--dirty'],
                                          cwd=work_dir,
                                          universal_newlines=True).strip()
    extra = ""
    if "dirty" in git_version:
        git_version = git_version.strip("-dirty")
        extra = "(working directory contains modifications)"
    print("Git repo commit id:\n\t", git_version, "\t", extra)

    # Add conda package versions


if __name__ == '__main__':
    # Create output directory from randomly generated name
    work_dir = os.path.dirname(os.path.abspath(__file__))

    # Call the click groups.
    main()

# snakemake -rp -s main.smk --cluster-config config/cluster.yaml --profile qsub_profile
# snakemake --dag -s main.smk| dot -Tpng > dag.png

# print dag functionality

