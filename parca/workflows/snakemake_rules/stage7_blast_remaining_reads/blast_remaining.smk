

checkpoint prepare_nt_blast_input:
    """
    Rule for creating blast multifasta infiles containing the number of reads from chunk_size.
    """
    input: 
        kmer_input="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage2/kmer_input/kmer_input.fasta",
        kmer_blast="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage6/kmer_species_subsetblast_classed.txt",
    output: 
        blast_infiles=directory("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage7/ntblastin"),
        read_count="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage7/count_printedfiles_assembledreadlengths.txt"
    params: 
        chunk_size=10000
    conda: "../../conda/biopython_env.yaml" #config['conda_environment'] 
    script:
        "../../scripts/blast_processing/blast_preprocessing/create_ntblast_input.py"

rule blast_remaining_reads_nt:
    """
    
    """
    input:
        infile="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage7/ntblastin/ntblast__{sliceiter}.fasta"
    output:
        blast_out="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage7/ntblastout/ntblast__{sliceiter}"
    params:
        blast_out_fmt="6 qseqid stitle evalue length nident sseqid bitscore staxids",
        e_val="1e-5",
        max_seqs="10",
        qcov_hsp_perc="20",
        nt_db_dir=config['nt_db_dir'] #config['nt_db_dir']
    conda: "../../conda/blast_env.yaml" #config['conda_environment'] 
    threads: 10
    shell:
        """
        blastn \
            -num_threads {threads} \
            -db {params.nt_db_dir}/nt \
            -query {input.infile} \
            -outfmt \'{params.blast_out_fmt}\' \
            -out {output.blast_out} \
            -evalue {params.e_val} \
            -max_target_seqs {params.max_seqs} \
            -qcov_hsp_perc {params.qcov_hsp_perc};
        """


def aggregate_ntblast_input(wildcards):
    checkpoint_output = checkpoints.prepare_nt_blast_input.get(**wildcards).output['blast_infiles']
    SLICEITER=glob_wildcards(os.path.join(checkpoint_output, "ntblast__{sliceiter, \d+}.fasta")).sliceiter

    slice_blast_list=expand("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage7/ntblastout/ntblast__{sliceiter}", 
           outdir=wildcards.outdir,
           sample=wildcards.sample,
           sample_type=wildcards.sample_type,
           nucleotide=wildcards.nucleotide,
           sliceiter=SLICEITER)

    return slice_blast_list

rule merge_nt_blast_result:
    """
    Rule for merging all blasted slices into one file, where a reads best blast result is kept.
    Input:
        blast_output = All blast result files from the nt_blast.
    Output:
        best_blast = Filtered blast results file, with the best blast hit.
            (Column names: "qseqid", "staxids", "score")
    """
    input: 
        blast_output=aggregate_ntblast_input
    output: 
        best_blast="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage7/best_nt_blast.txt"
    conda: "../../conda/R_env.yaml" #config['conda_environment'] 
    script:
        "../../scripts/blast_processing/blast_postprocessing/merge_slice_blast_result.R"

rule taxonomic_lineage_best_nt_blast:
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
        best_blast="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage7/best_nt_blast.txt"
    output:
        tax_id_lineage="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage7/best_nt_blast_tax_id_lineage.txt"
    conda: "../../conda/taxonkit_env.yaml" #config['conda_environment']
    params:
        dmp_dir=config['names_nodes_dmp_dir'] #config['names_nodes_dmp_dir']
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

rule reformat_nt_blast_taxids:
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
        best_blast="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage7/best_nt_blast.txt",
        tax_id_lineage="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage7/best_nt_blast_tax_id_lineage.txt"
    output: 
        species_and_above="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage7/best_nt_blast_species_and_above.txt",
        count_reads_tax_ids="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage7/count_reads_taxid_BLASTnt.txt"
    conda: "../../conda/R_env.yaml" #config['conda_environment'] 
    params:
        blast_type="BLASTnt",
        primates_file=config['primates_file']
    script:
        "../../scripts/blast_processing/blast_postprocessing/add_lineage_species_and_above.R"

rule merge_kmer_blast_nt_classed:
    input: 
        kmer_blast="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage6/kmer_species_subsetblast_classed.txt",
        species_and_above="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage7/best_nt_blast_species_and_above.txt"
    output: 
        kmer_blast_nt="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage7/kmer_species_subsetblast_blastnt_classed.txt",
        count_kmer_blast="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage7/count_kmer_SubsetBLAST_BLASTnt.txt"
    shell: 
        """
        cat {input.kmer_blast}  > {output.kmer_blast_nt};
        awk '$1!="seq_id"' {input.species_and_above} >> {output.kmer_blast_nt};
        count=$(awk '$1!="seq_id"' {output.kmer_blast_nt} | wc -l|cut -d " " -f 1);
        echo "type\tcount\nreads\t$count" > {output.count_kmer_blast};
        """
