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
    input:
        combined="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/combined_kraken_kaiju.txt"
    output:
        named="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/combined_kraken_kaiju_names_unfiltered.txt"
    params:
        names_nodes_dmp_dir=config['names_nodes_dmp_dir']
        #nodes=config['nodes'],
        #names=config['names']
    conda: config['conda_environment']
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
    input: 
        combined_unfiltered="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/combined_kraken_kaiju_names_unfiltered.txt",
        kraken_doublets="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/kraken_doublets.txt",
        kaiju_doublets="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/kaiju_doublets.txt",
        singletons="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/comparison/singletons.txt",
    output: 
        combined_SGF_empty_filter="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/combined_kraken_kaiju_names.txt",
        singletons_added_SGF_empty="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/singletons_added_SGF_empty.txt"
    conda: config['conda_environment']
    script:
        "../../scripts/taxonomy_processing/filter_SGF_empty.R" 

rule add_taxonomic_lineage_singletons:
    input: 
        singletons_added_SGF_empty="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/singletons_added_SGF_empty.txt",
    output:
        tax_id_lineage="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/singletons_tax_id_lineage.txt"
    conda: config['conda_environment']
    params:
        dmp_dir=config['names_nodes_dmp_dir']
    shell:
        """
        cut -f 3 {input.singletons_added_SGF_empty} | \
        awk '$1!="tax_id"' | sort | uniq |
        taxonkit lineage \
            --data-dir {params.dmp_dir} \
            --threads 2 --show-rank | \
        taxonkit reformat \
            --data-dir {params.dmp_dir} \
            --show-lineage-taxids | \
        cut -f 1,3,5 > {output.tax_id_lineage}
        """ 

rule singletons_species_to_genus:
    input: 
        singletons_added_SGF_empty="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/singletons_added_SGF_empty.txt",
        tax_id_lineage="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/singletons_tax_id_lineage.txt"
    output: 
        singletons_genus="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/singletons_genus.txt",
        read_count_singletons_genus = "{outdir}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage4/count_singletons_genus.txt"
    conda: config['conda_environment']
    script: 
        "../../scripts/taxonomy_processing/filter_taxonomy.R" 

rule add_taxon_names_singletons:
    input: 
        singletons_genus="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/singletons_genus.txt",
    output: 
        singletons_genus_names="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/singletons_genus_names.txt",
    params:
        names_nodes_dmp_dir=config['names_nodes_dmp_dir']
        #nodes=config['nodes'],
        #names=config['names']
    conda: config['conda_environment']
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
    input: 
        singletons_genus_names="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/singletons_genus_names.txt"
    output: 
        singletons_genus_names_reformat="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/singletons_genus_names_reformat.txt"
    conda: config['conda_environment']
    script:
        "../../scripts/taxonomy_processing/reformat_taxonomy.R"  

rule merge_combined_with_singletons:
    input: 
        singletons_genus_names_reformat="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/singletons_genus_names_reformat.txt",
        combined_SGF_empty_filter="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/combined_kraken_kaiju_names.txt"
    output: 
        combined_doublets_singletons="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/taxonomy_processing/combined_doublets_singletons.txt"
    shell: 
        """
        cat {input.singletons_genus_names_reformat} > {output.combined_doublets_singletons}
        awk '$1!="classified"' {input.combined_SGF_empty_filter} >> {output.combined_doublets_singletons}
        """