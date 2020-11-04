# rule create_taxonomy_db:
#     input:
#         names_dmp=config['names'],
#         nodes_dmp=config['nodes']
#     output:
#         sql_db="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/fiter_combined/accessionTaxa.sql"
#     conda: config['conda_environment']
#     script:
#         "../../scripts/taxonomy_processing/create_taxonomy_db.R"
rule add_taxon_names_doublets:
    """ 
    Rule for adding lineage to taxonomic IDs.
    Input: 
        combined=The merged kraken and kaiju results for the sequences that both softwares could classify.
    Params: 
        names_nodes_dmp_dir=config['names_nodes_dmp_dir']
    Output: 
        named=The lineage added to the input file.
    """ 
    input:
        combined="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/combined_kraken_kaiju.txt"
    output:
        named="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/combined_kraken_kaiju_names_unfiltered.txt"
    params:
        names_nodes_dmp_dir=config['names_nodes_dmp_dir'] #config['names_nodes_dmp_dir']
        #nodes=config['nodes'],
        #names=config['names']
    conda: "../../conda/kaiju_env.yaml" #config['conda_environment']
    shell:
        """
        kaiju-addTaxonNames \
            -t {params.names_nodes_dmp_dir}/nodes.dmp \
            -n {params.names_nodes_dmp_dir}/names.dmp \
            -i {input.combined} \
            -o {output.named} \
            -r superkingdom,class,order,family,genus,species;
        """


rule filter_SGF_empty:
    """ 
    Rule for filtering doublets sequence classifications (i.e. sequences that could be classified by both kraken and kaiju) to only contain classifications were species AND genus AND family is not NA.
    Doublets that are NA for species AND genus AND family are added to singletons. 
    Doublets that had either species OR genus OR family was not added to a file.

    Input: 
        combined_unfiltered=The merged kraken and kaiju results for the sequences that both softwares could classify where the lineage is added.
        kraken_doublets=Kraken classifications that were also classified by Kaiju.
        kaiju_doublets=Kaiju classifications that were also classified by Kraken.
        singletons=Classifications that could be made by only Kraken or Kaiju.
    Output: 
        combined_SGF_empty_filter=Classifications where species AND genus AND family exists.
        singletons_added_SGF_empty=Classifications made by kraken or kaiju and classifications that was missing species AND genus AND family.
    """ 
    input: 
        combined_unfiltered="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/combined_kraken_kaiju_names_unfiltered.txt",
        kraken_doublets="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/kraken_doublets.txt",
        kaiju_doublets="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/kaiju_doublets.txt",
        singletons="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/singletons.txt",
    output: 
        combined_SGF_empty_filter="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/combined_kraken_kaiju_names.txt",
        singletons_added_SGF_empty="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/singletons_added_SGF_empty.txt"
    conda: "../../conda/R_env.yaml" #config['conda_environment']
    script:
        "../../scripts/taxonomy_processing/filter_SGF_empty.R" 

rule get_taxonomic_lineage_singletons:
    """ 
    Rule for creating a dataframe with the taxonomic lineage for all taxids in the singletons file.
    Input: 
        singletons_added_SGF_empty=Classifications made by kraken or kaiju and classifications that was missing species AND genus AND family.
    Params: 
        dmp_dir=directory with names and nodes file
    Output: 
        tax_id_lineage=A dataframe with the taxonomic lineage for all taxids in singletons_added_SGF_empty.
    """ 
    input: 
        singletons_added_SGF_empty="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/singletons_added_SGF_empty.txt",
    output:
        tax_id_lineage="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/singletons_tax_id_lineage.txt"
    conda: "../../conda/taxonkit_env.yaml" #config['conda_environment']
    params:
        dmp_dir=config['names_nodes_dmp_dir'] #config['names_nodes_dmp_dir']
    shell:
        """
        [ ! -s {input.singletons_added_SGF_empty} ] && touch {output.tax_id_lineage} || \
        cut -f 3 {input.singletons_added_SGF_empty} | \
        awk '$1!="tax_id"' | sort | uniq | \
        taxonkit lineage \
            --data-dir {params.dmp_dir} \
            --show-rank | \
        taxonkit reformat \
            --data-dir {params.dmp_dir} \
            --show-lineage-taxids | \
        cut -f 1,3,5 > {output.tax_id_lineage};
        """ 

rule singletons_species_to_genus:
    """ 
    Rule for setting ranks below genus to genus.
    Input: 
        singletons_added_SGF_empty=Classifications made by kraken or kaiju and classifications that was missing species AND genus AND family.
        tax_id_lineage=A dataframe with the taxonomic lineage for all taxids in singletons_added_SGF_empty.
    Output: 
        singletons_genus=classifications where ranks below genus are moved to genus.
    """ 
    input: 
        singletons_added_SGF_empty="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/singletons_added_SGF_empty.txt",
        tax_id_lineage="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/singletons_tax_id_lineage.txt"
    output: 
        singletons_genus="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/singletons_genus.txt",
        read_count_singletons_genus = "{outdir}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage4/count_singletons_genus.txt"
    conda: "../../conda/R_env.yaml" #config['conda_environment']
    script: 
        "../../scripts/taxonomy_processing/filter_taxonomy.R" 

rule add_taxon_names_singletons:
    """ 
    Rule for 
    Input: 
    Params: 
    Output: 
    """ 
    input: 
        singletons_genus="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/singletons_genus.txt",
    output: 
        singletons_genus_names="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/singletons_genus_names.txt",
    params:
        names_nodes_dmp_dir=config['names_nodes_dmp_dir'] #config['names_nodes_dmp_dir']
        #nodes=config['nodes'],
        #names=config['names']
    conda: "../../conda/kaiju_env.yaml" #config['conda_environment']
    shell:
        """
        kaiju-addTaxonNames \
        -t {params.names_nodes_dmp_dir}/nodes.dmp \
        -n {params.names_nodes_dmp_dir}/names.dmp \
        -i {input.singletons_genus} \
        -o {output.singletons_genus_names} \
        -r superkingdom,class,order,family,genus,species;
        """

rule reformat_singletons:
    """ 
    Rule for 
    Input: 
    Params: 
    Output: 
    """ 
    input: 
        singletons_genus_names="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/singletons_genus_names.txt"
    output: 
        singletons_genus_names_reformat="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/singletons_genus_names_reformat.txt"
    conda: "../../conda/R_env.yaml" #config['conda_environment']
    script:
        "../../scripts/taxonomy_processing/reformat_taxonomy.R"  

rule merge_combined_with_singletons:
    """ 
    Rule for 
    Input: 
    Params: 
    Output: 
    """ 
    input: 
        singletons_genus_names_reformat="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/singletons_genus_names_reformat.txt",
        combined_SGF_empty_filter="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/combined_kraken_kaiju_names.txt"
    output: 
        combined_doublets_singletons="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/combined_doublets_singletons.txt"
    shell: 
        """
        cat {input.singletons_genus_names_reformat} > {output.combined_doublets_singletons};
        awk '$1!="classified"' {input.combined_SGF_empty_filter} >> {output.combined_doublets_singletons};
        """

rule genus_species_split:
    """ 
    Rule for 
    Input: 
    Params: 
    Output: 
    """ 
    input: 
        combined_doublets_singletons="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/combined_doublets_singletons.txt"
    output:
        species="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/genusspeciessplit/species_classed.txt",
        higher="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/genusspeciessplit/above_species_classed.txt",
        read_count = "{outdir}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage4/count_species_genus_higher.txt"
    conda: "../../conda/R_env.yaml" #config['conda_environment']
    script: "../../scripts/taxonomy_processing/genus_species_split.R" 