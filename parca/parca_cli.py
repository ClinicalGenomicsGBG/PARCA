
import os
import click
import shutil
import snakemake
import subprocess
import pandas as pd
from workflows.utils.process_runinfo_metadata import ProcessRuninfoMetadata
from snakemake.logging import logger


@click.group()
def main():
    """
    P a R C A - Pathogen Research in Clinical Applications
    """
    pass


@main.command()
@click.option('-m', '--metadata', 'metadata',
              required=True,
              type=click.Path(exists=True))
@click.option('-r', '--runinfo', 'runinfo',
              required=True,
              type=click.Path(exists=True))
@click.option('-o', '--outdir', 'outdir', type=click.Path(exists=True),
              required=True,
              help='Give the absolute path to a directory where all results will be placed in a subforder named after date and runinfo')
@click.option('-cl', '--complete_log', 'complete_log',
              type=click.Path(exists=True),
              default="/medstore/logs/pipeline_logfiles/parca",
              help='Path to a log directory')
# @click.option('-gs', '--generate_subdir', 'generate_subdir', is_flag=True,
#               help='Generate a subfolder with date within outdir a given outdir')
@click.option('-d', '--dryrun', 'dryrun', is_flag=True, help='dryrun')
def run(metadata, runinfo, dryrun, outdir, complete_log):
    """
    Run the PaRCA pipeline.
    """
    # if generate_subdir and outdir:
    #     base_outdir = outdir
    #     sub_outdir = GenerateOutdir.get_date_and_randomizer()
    #     outdir = os.path.join(base_outdir, sub_outdir)

    run_dict_list = ProcessRuninfoMetadata.generate_runinfo_dict(runinfo)
    metadata_dict = ProcessRuninfoMetadata.generate_metadata_dict(metadata)

    config_dict_added = {
                'run_dict_list': run_dict_list,
                'metadata_dict': metadata_dict,
                'outdir': outdir}

    cluster_settings = "".join(["qsub ",
                                "-S /bin/bash ",
                                "-pe mpi {cluster.threads} ",
                                "-q {cluster.queue} ",
                                "-S /bin/bash ",
                                "-N parca-{rule}-{wildcards.run_id} ",
                                "-V ",
                                "-cwd ",
                                "-j y -o {log}.cluster"])  #-l excl=1

    status = snakemake.snakemake(snakefile=f'{work_dir}/main.smk',
                                 dryrun=dryrun,
                                 config=config_dict_added,
                                 workdir=work_dir,
                                 latency_wait=60,
                                 printreason=True,
                                 printshellcmds=True,
                                 verbose=True,
                                 cores=40,
                                 # conda settings
                                 use_conda=True,
                                 conda_prefix=f'{outdir}/conda',
                                 # conda_cleanup_envs=True,
                                 # conda_create_envs_only=True,
                                 # Singularity settings # singularity_args=" --cleanenv "
                                 use_singularity=True)
                                 # cluster settings
                                 #cluster_config=f'{work_dir}/config/cluster.yaml',
                                 #cluster=cluster_settings,
                                 #max_jobs_per_second=99)

                                 # Double check if this can be replaced with qsub profile... could not find this...
                                 #  cleanup_shadow=True)
                                 #  conda_cleanup_envs=True)

    print("STATUSCODE:", status)  # True or False

    if not dryrun:
        run_ids = "_".join([f"{run_dict.get('start_date')}-{run_dict.get('run_id')}" for run_dict in run_dict_list])

        logfile_path = logger.logfile
        new_log_dst = os.path.join(complete_log, f'{run_ids}_snakemake.log')
        if os.path.exists(new_log_dst):
            os.remove(new_log_dst)

        shutil.copy(logfile_path, new_log_dst)


    # clean up if error... do not use this since if the generate outdir is not used it will remove unnecessary things 
    # return_code=subprocess.call(['rmdir', outdir])


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


# print dag functionality
# snakemake --dag -s main.smk| dot -Tpng > dag.png

# python3 parca_cli.py run -m /apps/bio/dev_repos/parca/demo/runinfo/metadata.csv -r /apps/bio/dev_repos/parca/demo/runinfo/runinfo.csv -o /medstore/logs/pipeline_logfiles/parca --dryrun