
checkpoint tableview:
    input: 
        read_count="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/readcount.tsv",
        trimmed_read_count="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage1/trimming/count_bbduk_trimmed_reads.txt"
    output: 
        tableview="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/readcount_tableview.tsv",
        classified_reads_mincount="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/classified_reads_mincount.tsv",
        read_count="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage8/tableview/count_classified_reads.txt",
        organism_dir=directory("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/organism_dir"),
        kingdom_dir=directory("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/kingdom_dir")
    params:
        #SE_or_PE="SE",
        mincount=config['tableview_min_count']
    conda: "../../conda/R_env.yaml"
    script: "../../scripts/reformat_results/tableview_splitting.R"

rule filter_fastq_organism_SE:
    input: 
        organism_tableview="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage8/tableview/organism_dir/taxid_{taxid}.tsv",
        reads="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage1/trimming/trimmed_reads.fq"
    output: 
        fastq_out="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage8/tableview/organism_fastq/{taxid}.fastq"
    params:
        SE_or_PE="SE",
        negate_query="FALSE"  # The string should be in uppercase letters for R to interpret this.
    conda: "../../conda/R_env.yaml"
    script: "../../scripts/reformat_results/filter_fastq.R"

rule filter_fastq_organism_PE:
    input: 
        organism_tableview="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage8/tableview/organism_fastq/taxid_{taxid}.tsv",
        unmerged_reads="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/unmerged_reads_trimmed.fq",
        merged_reads="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/merged_reads_trimmed.fq"
    output: 
        fastq_out="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage8/tableview/organism_fastq/{taxid}.fastq"
    params:
        SE_or_PE="PE",
        negate_query="FALSE"  # The string should be in uppercase letters for R to interpret this.
    conda: "../../conda/R_env.yaml"
    script: "../../scripts/reformat_results/filter_fastq.R"


rule filter_fastq_kingdom_SE:
    input: 
        organism_tableview="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage8/tableview/kingdom_dir/kingdom_{kingdom}.tsv",
        reads="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage1/trimming/trimmed_reads.fq"
    output: 
        fastq_out="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage8/tableview/kingdom_fastq/{kingdom}.fastq"
    params:
        SE_or_PE="SE",
        negate_query="FALSE"  # The string should be in uppercase letters for R to interpret this.
    conda: "../../conda/R_env.yaml"
    script: "../../scripts/reformat_results/filter_fastq.R"

rule filter_fastq_kingdom_PE:
    input: 
        organism_tableview="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage8/tableview/kingdom_fastq/kingdom_{kingdom}.tsv",
        unmerged_reads="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/unmerged_reads_trimmed.fq",
        merged_reads="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/merged_reads_trimmed.fq"
    output: 
        fastq_out="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage8/tableview/kingdom_fastq/{kingdom}.fastq"
    params:
        SE_or_PE="PE",
        negate_query="FALSE"  # The string should be in uppercase letters for R to interpret this.
    conda: "../../conda/R_env.yaml"
    script: "../../scripts/reformat_results/filter_fastq.R"

rule filter_fastq_unclassified_SE:
    input: 
        organism_tableview="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage8/tableview/classified_reads_mincount.tsv",
        reads="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage1/trimming/trimmed_reads.fq"
    output: 
        fastq_out="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage8/tableview/unclassified_fastq/unclassified.fastq"
    params:
        SE_or_PE="SE",
        negate_query="TRUE"  # The string should be in uppercase letters for R to interpret this.
    conda: "../../conda/R_env.yaml"
    script: "../../scripts/reformat_results/filter_fastq.R"

rule filter_fastq_unclassified_PE:
    input: 
        organism_tableview="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage8/tableview/classified_reads_mincount.tsv",
        unmerged_reads="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/unmerged_reads_trimmed.fq",
        merged_reads="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage1/trimming/merged_reads_trimmed.fq"
    output: 
        fastq_out="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/PE_{nucleotide}/stage8/tableview/unclassified_fastq/unclassified.fastq"
    params:
        SE_or_PE="PE",
        negate_query="TRUE"  # The string should be in uppercase letters for R to interpret this.
    conda: "../../conda/R_env.yaml"
    script: "../../scripts/reformat_results/filter_fastq.R"

rule zip_filtered_fastq_organism:
    input: 
        fastq="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/organism_fastq/{taxid}.fastq"
    output: 
        fastq="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/organism_fastq/{taxid}.fastq.gz"
    threads: 4
    shell:
        """
        pigz -p {threads} -k {input.fastq};
        """ 

rule zip_filtered_fastq_kingdom:
    input: 
        fastq="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/kingdom_fastq/{kingdom}.fastq"
    output: 
        fastq="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/kingdom_fastq/{kingdom}.fastq.gz"
    threads: 4
    shell:
        """
        pigz -p {threads} -k {input.fastq};
        """ 

rule zip_filtered_fastq_unclassified:
    input: 
        fastq="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}/stage8/tableview/unclassified_fastq/unclassified.fastq"
    output: 
        fastq="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}stage8/tableview/unclassified_fastq/unclassified.fastq.gz",
        read_count="{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/stats_{sample_type}_{nucleotide}/stage8/tableview/count_unclassified_reads.txt"
    threads: 8
    shell:
        """
        echo count > {output.read_count};
        echo $(cat {input.fastq}|wc -l)/4|bc  >> {output.read_count};
        pigz -p {threads} -k {input.fastq};
        """ 

def filter_fastq_according_to_classification(wildcards):
    checkpoint_output_organism = checkpoints.tableview.get(**wildcards).output['organism_dir']
    organism_list = expand("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/stage8/tableview/organism_fastq/{taxid}.fastq.gz",
           outdir=wildcards.outdir,
           start_date=wildcards.start_date,
           run_id=wildcards.run_id,
           sample=wildcards.sample,
           nucleotide=wildcards.nucleotide,
           taxid=glob_wildcards(os.path.join(checkpoint_output_organism, "{taxid, \d+}")).taxid)

    checkpoint_output_kingdom = checkpoints.tableview_SE.get(**wildcards).output['kingdom_dir']
    kingdom_list = expand("{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/SE_{nucleotide}/tableview/kingdom_fastq/{kingdom}.fastq.gz",
           outdir=wildcards.outdir,
           start_date=wildcards.start_date,
           run_id=wildcards.run_id,
           sample=wildcards.sample,
           nucleotide=wildcards.nucleotide,
           kingdom=glob_wildcards(os.path.join(checkpoint_output_kingdom, "{kingdom, \d+}")).kingdom)

    organism_kingdom_list = organism_list + kingdom_list
    return organism_kingdom_list

rule call_filter_fastqs:
    input:
        filter_fastq_according_to_classification,
        "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}stage8/tableview/unclassified_fastq/unclassified.fastq.gz"
    output:
        "{outdir}/{start_date}_{run_id}/snakemake_results_{sample}/{sample_type}_{nucleotide}stage8/tableview/fastq_filtering_done"
    shell:
        "touch {output}"

