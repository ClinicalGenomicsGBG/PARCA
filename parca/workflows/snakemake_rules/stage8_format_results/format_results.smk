rule rename_classified:
    input: 
        kmer_blast_nt="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage7/kmer_species_subsetblast_blastnt_classed.txt",
    output: 
        all_classed="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/all_classed_scores.txt"
    shell: 
        """
        cp {input.kmer_blast_nt} {output.all_classed}
        """

rule format_all_classified:
    input: 
        all_classed="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/all_classed_scores.txt"
    output: 
        all_classed_read_taxid="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/all_classed_read_taxid.txt",
        type_summary="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/type_summary.txt"
    conda: "../../conda/R_env.yaml" #config['conda_environment'] 
    script:
         "../../scripts/reformat_results/reformat_all_classed.R"

rule add_taxon_names_all_classed:
    input: 
        all_classed_read_taxid="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/all_classed_read_taxid.txt"
    output: 
        all_classed_read_taxid_names="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/all_classed_read_taxid_names.txt"
    conda: "../../conda/kaiju_env.yaml" #config['conda_environment'] 
    params: 
        names_nodes_dmp_dir=config['names_nodes_dmp_dir'] #config['names_nodes_dmp_dir']
    shell:
        """
        kaiju-addTaxonNames \
        -t {params.names_nodes_dmp_dir}/nodes.dmp \
        -n {params.names_nodes_dmp_dir}/names.dmp \
        -i {input.all_classed_read_taxid} \
        -o {output.all_classed_read_taxid_names} \
        -r superkingdom,phylum,order,family,genus,species;
        """ 

#samtools view -b reads.bam chr1:10420000-10421000 > subset.bam