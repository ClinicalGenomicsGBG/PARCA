
rule existing_slices_split:
    input: 
        higher="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/genusspeciessplit/above_species_classed.txt",
    output: 
        detected="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/downloadblastslices/detected_slices.txt",
        missing="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/downloadblastslices/missing_slices.txt",
        read_count="{outdir}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage5/count_detected_missing_slices.txt"
    params: 
        existing_slice_path=runinfo_dict['existing_slice_path'], #config['existing_slice_path'],
        min_tax_id_count=2
    conda: "../../conda/R_env.yaml" #config['conda_environment']
    script: "../../scripts/blast_processing/blast_preprocessing/existing_slices_split.R" 

def create_file_list(file_name):
    with open(file_name, 'r') as slice_file:
        slice_name = [line.strip() for line in slice_file.readlines()]
    return(slice_name)

rule copy_detected_slices:
    input: 
        detected="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/downloadblastslices/detected_slices.txt",
    output: temp(directory("{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/downloadblastslices/existing_slices"))
    params: 
        existing_slice_path=runinfo_dict['existing_slice_path'] #config['existing_slice_path']
    benchmark: "{outdir}/snakemake_results_{sample}/benchmarks_{sample_type}_{nucleotide}/stage5/existing_slices_copy.txt"
    run: 
        detected_list=create_file_list(input.detected)
        shell("if [ ! -d {output} ]; then \mkdir {output};fi;")
        for slice_file in detected_list:
            shell("cp {params.existing_slice_path}/{slice_file} {output}/{slice_file}")


rule download_missing_slices:
    input: 
        missing="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/downloadblastslices/missing_slices.txt"
    output: temp(directory("{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/downloadblastslices/downloaded_slices"))
    params: 
        existing_slice_path=runinfo_dict['existing_slice_path'] #config['existing_slice_path']
    benchmark: "{outdir}/snakemake_results_{sample}/benchmarks_{sample_type}_{nucleotide}/stage5/download_slices.txt"
    run: 
        missing_list=create_file_list(input.missing)
        shell("if [ ! -d {output} ]; then \mkdir {output};fi;")
        if len(missing_list) > 0:
            for slice_file in missing_list:
                shell("wget -q -t 20 --waitretry=10 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=txid{slice_file}[Subtree]&retmax=400000&rettype=uilist' -O {output}/{slice_file}")

rule all_downloaded_slices:
    input: 
        downloaded_dir="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/downloadblastslices/downloaded_slices"
    output: 
        all_downloaded=temp("{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/downloadblastslices/downloaded_slices.txt")
    conda: "../../conda/R_env.yaml" #config['conda_environment']
    benchmark: "{outdir}/snakemake_results_{sample}/benchmarks_{sample_type}_{nucleotide}/stage5/download_slices_merged.txt"
    script: "../../scripts/blast_processing/blast_preprocessing/merge_downloaded_slices.R" 

rule create_tax_id_accession_slice_files:
    """
    Rule for filtering the downloaded slices from a dmp file.
    Returns a directory with files called a after a taxid containing all accessions for that taxid.
    """
    input: 
        all_downloaded="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/downloadblastslices/downloaded_slices.txt"
    output: 
        created_slice_dir=temp(directory("{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/downloadblastslices/downloaded_slices_acc"))
    params:
        splitaccdump_dir=runinfo_dict['splitaccdump_dir'] #config['splitaccdump_dir']
    conda: "../../conda/R_env.yaml" #config['conda_environment'] 
    benchmark: "{outdir}/snakemake_results_{sample}/benchmarks_{sample_type}_{nucleotide}/stage5/downloaded_slices_acc.txt"
    script:  "../../scripts/blast_processing/blast_preprocessing/create_slice_files_downloaded.R"

checkpoint all_gislices:
    input: 
        created_slice_dir="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/downloadblastslices/downloaded_slices_acc",
        existing_slices="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/downloadblastslices/existing_slices"
    output: 
        all_slices=directory("{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/downloadblastslices/all_gislices")
    params: 
        splitaccdump_dir=runinfo_dict['splitaccdump_dir'] #config['splitaccdump_dir']
    shell: 
        """
        if [ ! -d {output.all_slices} ]; then \mkdir {output.all_slices};fi;

        if [ "$(ls {input.created_slice_dir})" ]; then
            cp {input.created_slice_dir}/* {output.all_slices};
            cp {input.created_slice_dir}/* {params.splitaccdump_dir};
        fi;

        if [ "$(ls {input.existing_slices})" ]; then
            cp {input.existing_slices}/* {output.all_slices};
        fi;
        """

rule create_blastdb_alias:
    input:
        "{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/downloadblastslices/all_gislices/{gi_slice}"
    output:
        "{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/blastslices/{gi_slice}.nal"
    params: 
        nt_db_dir=runinfo_dict['nt_db_dir'], #config['nt_db_dir'],
        out="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/blastslices/{gi_slice}"
    conda: "../../conda/blast_env.yaml" #config['conda_environment'] 
    shell:
        """
        blastdb_aliastool -dbtype nucl -seqidlist {input} -db {params.nt_db_dir}/nt -out {params.out} >/dev/null;
        """


def aggregate_blast_slices(wildcards):
    checkpoint_output = checkpoints.all_gislices.get(**wildcards).output[0]
    return expand("{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/blastslices/{gi_slice}.nal",
           outdir=wildcards.outdir,
           sample=wildcards.sample,
           sample_type=wildcards.sample_type,
           nucleotide=wildcards.nucleotide,
           gi_slice=glob_wildcards(os.path.join(checkpoint_output, "{gi_slice, \d+}")).gi_slice)


# an aggregation over all produced clusters
rule call_create_blastdb_alias:
    input:
        aggregate_blast_slices
    output:
        "{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/alias_done"
    shell:
        "touch {output}"


# checkpoint prepare_blast_input:
#     input: 
#         kmer_input="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage2/kmer_input/kmer_input.fasta",
#         higher="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/genusspeciessplit/above_species_classed.txt"
#     output: 
#         blast_infiles=directory("{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/sliceblastin"),
#         count="{outdir}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage5/count_printedfiles_assembledreadlengths.txt"
#     params: 
#         chunk_size=6000
#     conda: config['conda_environment'] 
#     script:
#         "../../scripts/blast_processing/blast_preprocessing/create_sliceblast_input.py"