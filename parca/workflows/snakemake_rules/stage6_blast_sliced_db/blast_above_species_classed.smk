
checkpoint prepare_blast_input:
    input: 
        kmer_input="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage2/kmer_input/kmer_input.fasta",
        higher="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/genusspeciessplit/above_species_classed.txt"
    output: 
        blast_infiles=directory("{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/sliceblastin"),
        count="{outdir}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage5/count_printedfiles_assembledreadlengths.txt"
    params: 
        chunk_size=6000
    conda: config['conda_environment'] 
    script:
        "../../scripts/blast_processing/blast_preprocessing/create_sliceblast_input.py"


rule blast_slices:
    input:
        infile_slice="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/sliceblastin/{gi_slice}__{sliceiter}",
        db_slice="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/blastslices/{gi_slice}.nal"
    output:
        blast_out="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage6/sliceblastout/{gi_slice}__{sliceiter}"
    params:
        blast_out_fmt="6 qseqid stitle evalue length nident sseqid bitscore staxids",
        e_val="1e-3",
        max_seqs="10",
        db_dir="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage5/blastslices/{gi_slice}"
    conda: config['conda_environment'] 
    threads: 10
    shell:
        """
        blastn \
            -num_threads {threads} \
            -db {params.db_dir} \
            -query {input.infile_slice} \
            -outfmt \'{params.blast_out_fmt}\' \
            -out {output.blast_out} \
            -evalue {params.e_val} \
            -max_target_seqs {params.max_seqs};
        """

def aggregate_sliceblast_input(wildcards):
    checkpoint_output = checkpoints.prepare_blast_input.get(**wildcards).output['blast_infiles']
    wildcard_obj=glob_wildcards(os.path.join(checkpoint_output, "{gi_slice, \d+}__{sliceiter, \d+}"))
    GI_SLICE=wildcard_obj.gi_slice
    SLICEITER=wildcard_obj.sliceiter

    slice_blast_list=expand("{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage6/sliceblastout/{gi_slice}__{sliceiter}", zip,
           outdir=[wildcards.outdir]*len(GI_SLICE),
           sample=[wildcards.sample]*len(GI_SLICE),
           sample_type=[wildcards.sample_type]*len(GI_SLICE),
           nucleotide=[wildcards.nucleotide]*len(GI_SLICE),
           gi_slice=GI_SLICE,
           sliceiter=SLICEITER)

    return slice_blast_list

rule merge_slice_blast_result:
    """
    Rule for merging all blasted slices into one file, where a reads best blast result is kept.
    Input:
        blast_output = All blast result files from the blasting of slices.
    Output:
        best_blast = Filtered blast results file, with the best blast hit.
            (Column names: "qseqid", "staxids", "score")
    """
    input: 
        blast_output=aggregate_sliceblast_input
    output: 
        best_blast="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage6/best_blast.txt"
    conda: config['conda_environment'] 
    script:
        "../../scripts/blast_processing/blast_postprocessing/merge_slice_blast_result.R"

rule taxonomic_lineage_best_blast:
    """
    Rule for adding the taxonomic lineage (above a current taxid) to the taxids from best blast
    Input:
        best_blast = Filtered blast results file, with the best blast hit.
            (Column names: "qseqid", "staxids", "score")
    Params:
        dmp_dir = Name of directory where the names.dmp and nodes.dmp are located.
    Output:
        tax_id_lineage = Lineage file of taxids above a current taxid.
            (Columns: taxid, rank, taxid lineage separated by ";")
    """
    input: 
        best_blast="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage6/best_blast.txt"
    output:
        tax_id_lineage="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage6/best_blast_tax_id_lineage.txt"
    conda: config['conda_environment']
    params:
        dmp_dir=config['names_nodes_dmp_dir']
    shell:
        """
        [ ! -s {input.best_blast} ] && touch {output.tax_id_lineage} || \
        cut -f 2 {input.best_blast} | \
        awk '$1!="staxids"' | sort | uniq | \
        taxonkit lineage \
            --data-dir {params.dmp_dir} \
            --show-rank | \
        taxonkit reformat \
            --data-dir {params.dmp_dir} \
            --show-lineage-taxids | \
        cut -f 1,3,5 > {output.tax_id_lineage};
        """ 

rule reformat_blast_taxids:
    """
    Rule for moving all reads classed to taxids below species, up to species.
    Input:
        best_blast = Filtered blast results file, with the best blast hit.
            (Column names: "qseqid", "staxids", "score")
        tax_id_lineage = Lineage file of taxids above a current taxid.
            (Columns: taxid, rank, taxid lineage separated by ";")
    Params:
        blast_type = Name of the analysis e.g. SubsetBLAST, the type column of species_and_above output.
        primates_file = xml file with taxids that should be classed as primates.
    Output:
        species_and_above = Outfile where reads classed below species are moved up to species.
            (Column names: "seq_id", "tax_id", "score", "rank", "organism", "type")
        count_reads_tax_ids = Stats of how many unique taxids are found and how many reads that are classed.
    """
    input: 
        best_blast="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage6/best_blast.txt",
        tax_id_lineage="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage6/best_blast_tax_id_lineage.txt"
    output: 
        species_and_above="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage6/best_blast_species_and_above.txt",
        count_reads_tax_ids="{outdir}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage6/count_reads_taxid_SubsetBLAST.txt"
    conda: config['conda_environment'] 
    params:
        blast_type="SubsetBLAST",
        primates_file=config['primates_file']
    script:
        "../../scripts/blast_processing/blast_postprocessing/add_lineage_species_and_above.R"

rule merge_kmer_classed_and_blast_classed:
    input: 
        species="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage4/genusspeciessplit/species_classed.txt",
        species_and_above="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage6/best_blast_species_and_above.txt"
    output: 
        kmer_blast="{outdir}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage6/kmer_species_subsetblast_classed.txt",
        count_kmer_blast="{outdir}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage6/count_kmer_SubsetBLAST.txt"
    shell: 
        """
        cat {input.species}  > {output.kmer_blast};
        awk '$1!="seq_id"' {input.species_and_above} >> {output.kmer_blast};
        count=$(awk '$1!="seq_id"' {output.kmer_blast} | wc -l|cut -d " " -f 1);
        echo "type\tcount\nreads\t$count" > {output.count_kmer_blast};
        """

