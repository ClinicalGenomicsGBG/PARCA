
rule existing_slices_split:
    input: 
        higher="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/genusspeciessplit/above_species_classed.txt",
    output: 
        detected="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/downloadblastslices/detected_slices.txt",
        missing="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/downloadblastslices/missing_slices.txt",
        count="{outdir}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage5/count_detected_missing_slices.txt"
    params: 
        existing_slice_path=config['existing_slice_path'],
        min_tax_id_count=2
    conda: config['conda_environment']
    script: "../../scripts/blast_processing/existing_slices_split.R" 

def create_file_list(file_name):
    with open(file_name, 'r') as slice_file:
        slice_name = [line.strip() for line in slice_file.readlines()]
    return(slice_name)

rule copy_detected_slices:
    input: 
        detected="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/downloadblastslices/detected_slices.txt",
    output: directory("{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/downloadblastslices/existing_slices")
    params: 
        existing_slice_path=config['existing_slice_path']
    run: 
        detected_list=create_file_list(input.detected)
        shell("if [ ! -d {output} ]; then \mkdir {output};fi;")
        for slice_file in detected_list:
            shell("cp {params.existing_slice_path}/{slice_file} {output}/{slice_file}")


rule download_missing_slices:
    input: 
        missing="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/downloadblastslices/missing_slices.txt"
    output: directory("{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/downloadblastslices/downloaded_slices")
    params: 
        existing_slice_path=config['existing_slice_path']
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
        all_downloaded="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/downloadblastslices/downloaded_slices.txt"
     script: "../../scripts/blast_processing/existing_slices_split.R" 