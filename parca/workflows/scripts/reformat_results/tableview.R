# Author: Pernilla Ericsson (pernilla.ericsson@gu.se, clinicalgenomics@gu.se)
# Date: 2020-12-22

suppressPackageStartupMessages({
  library(tidyverse);library(magrittr);library(data.table);library(ShortRead)
} )

read_count <- "/medstore/logs/pipeline_logfiles/parca/20201202_run_1/snakemake_results_sample_1/PE_RNA/stage8/readcount.tsv"
fastq_path_merged <- "/medstore/logs/pipeline_logfiles/parca/20201202_run_1/snakemake_results_sample_1/PE_RNA/stage1/trimming/unmerged_reads_trimmed.fq"
fastq_path_unmerged <- "/medstore/logs/pipeline_logfiles/parca/20201202_run_1/snakemake_results_sample_1/PE_RNA/stage1/trimming/merged_reads_trimmed.fq"

trimmed_read_count_unmerged <-"/medstore/logs/pipeline_logfiles/parca/20201202_run_1/snakemake_results_sample_1/stats_PE_RNA/stage1/trimming/count_bbduk_unmerged_reads_trimmed_raw.txt"
trimmed_read_count_merged <- "/medstore/logs/pipeline_logfiles/parca/20201202_run_1/snakemake_results_sample_1/stats_PE_RNA/stage1/trimming/count_bbduk_merged_reads_trimmed.txt"

SE_or_PE <- "PE"
mincount <- 4

tableview_out <- "/home/xerpey/projects/parca_tests/tableview.tsv"
read_classifications_out <- "/home/xerpey/projects/parca_tests/classifications.tsv"
unclassified_reads_out <- "/home/xerpey/projects/parca_tests/unclassified.fastq"
organism_dir <- "/home/xerpey/projects/parca_tests/organism_dir"
kingdom_dir <- "/home/xerpey/projects/parca_tests/kingdom_dir"


# read_count <- snakemake@input[['read_count']]
# 
# SE_or_PE <- snakemake@params[['SE_or_PE']]
# mincount <- snakemake@params[['mincount']]
# 
# tableview_out <- snakemake@output[['tableview']]
# read_classifications_out <- snakemake@output[['read_classifications']]
# unclassified_reads_out <- snakemake@output[['unclassified_reads']]
# 
# organism_dir <- snakemake@output[['organism_dir']]
# kingdom_dir <- snakemake@output[['kingdom_dir']]

# Read in the fastq files and the read count stats with regard to PE or SE
if (SE_or_PE == "SE") {
  trimmed_read_count <- snakemake@input[['trimmed_read_count']]
  fastq_path <- snakemake@input[['reads']]
  
  fastqs=c(fastq_path)
  
  read_total <- read_tsv(trimmed_read_count) %>% pull(count)

} else if (SE_or_PE == "PE"){
  # trimmed_read_count_unmerged <- snakemake@input[['trimmed_read_count_unmerged']]
  # trimmed_read_count_merged <- snakemake@input[['trimmed_read_count_merged']]
  # 
  # fastq_path_merged <- snakemake@input[['unmerged_reads']]
  # fastq_path_unmerged <- snakemake@input[['merged_reads']]
  
  fastqs=c(fastq_path_merged,fastq_path_unmerged)
  
  read_total <- 
    c(trimmed_read_count_unmerged, trimmed_read_count_merged) %>% 
    map_dfr(~read_tsv(.x)) %>% select(count) %>% colSums() %>% unname()

} else {
  quit(save="no", status=0)
}

# Read in the classified reads dataframe with readcount included.
read_count_df <-data.table::fread(file = read_count, fill=TRUE, sep = "\t",
                                  header=TRUE) %>% as_tibble()

if(nrow(read_count_df)==0){
  if(!dir.exists(organism_dir)) dir.create(organism_dir, recursive = TRUE)
  if(!dir.exists(kingdom_dir)) dir.create(kingdom_dir, recursive = TRUE)
  write_tsv(tibble(), tableview_out)
  quit(save = "no", status = 0)
}

# Add a column with percent of total reads.
tableview <- 
  read_count_df %>% 
  group_by(superkingdom, organism, tax_id) %>% 
  summarize(read_count=sum(read_count, na.rm = TRUE), .groups = "drop")  %>% 
  mutate(percent=(read_count/read_total)*100, .before=read_count) %>% 
  filter(read_count>mincount) 

tableview %>% 
  write_tsv(tableview_out, col_names = TRUE)

# Write the reads ids with number of reads higher than mincount to a file.
unique_tax_id <- tableview %>% pull(tax_id)
organism_df_list <-
  read_count_df %>% select(tax_id, seq_id, superkingdom, organism) %>% 
  filter(tax_id %in% unique_tax_id) %>% write_tsv(read_classifications_out, col_names = TRUE)

append_to_fastq <- function(fastq, ids, outputFile) {
  # Function to append reads to a fastq file.
  i_query <- ids %>% map_chr(~{glue::glue("^{.x}\\s")}) %>% paste(collapse = "|")
  subset <- fastq[ str_detect(ShortRead::id(fastq), i_query ) ]
  print(outputFile)
  ShortRead::writeFastq(object = subset, file = outputFile, mode="a", compress = FALSE)
}

generate_subsetted_fastqs <- function(fastq_path, all_classed, unclassified_reads, outdir_organism, outdir_kingdom ) {
  # Function to call append_to_fastq and create a fastq file with unassigned reads.
  if(!dir.exists(outdir_organism)) dir.create(outdir_organism, recursive = TRUE)
  if(!dir.exists(outdir_kingdom)) dir.create(outdir_kingdom, recursive = TRUE)
  print("Outdirs created..")
  
  fastq <- ShortRead::readFastq(fastq_path)
  
  all_seq_id <- all_classed %>% pull(seq_id)
  all_query <- all_seq_id %>% map_chr(~{glue::glue("^{.x}\\s")}) %>% paste(collapse = "|")
  assigned_fastq <- fastq[ str_detect(ShortRead::id(fastq), all_query ) ]
  
  print("Working...")
  
  df_splitted <- all_classed %>% select(seq_id, tax_id) %$% 
    split(., tax_id) %>% 
    map( ~{pull(.x, seq_id)})
  print("Number of taxids:")
  print(length(df_splitted))
  imap(df_splitted, ~{
    print(glue::glue("Number of reads: {nrow(.x)}, Tax id: {.y} \n"))
    outfile_organism <- glue::glue("{outdir_organism}/taxid_{.y}.fastq")
    
    append_to_fastq(assigned_fastq, ids=.x, outputFile=outfile_organism)
  })
  
  df_splitted_kingdom <- all_classed %>% 
    mutate(superkingdom=str_trim(superkingdom) %>% str_replace_all(" ", "_") ) %>% 
    mutate(superkingdom=ifelse(is.na(superkingdom), "NA", superkingdom)) %>% 
    select(seq_id, superkingdom) %$% 
    split(., superkingdom) %>% 
    map( ~{pull(.x, seq_id)})
  
  print("Number if kingdoms")
  print(length(df_splitted_kingdom))
  
  imap(df_splitted_kingdom, ~{
    print(glue::glue("Number of reads: {nrow(.x)}, Kingdom: {.y} \n"))
    outfile_kingdom <- glue::glue("{outdir_kingdom}/kingdom_{.y}.fastq")
    
    append_to_fastq(assigned_fastq, ids=.x, outputFile=outfile_kingdom)
  })
  
  unassigned_fastq <- fastq[ str_detect(ShortRead::id(fastq), all_query, negate = TRUE) ]
  ShortRead::writeFastq(object = unassigned_fastq, file = unclassified_reads, mode="a", compress = TRUE)
}

all_classed <- organism_df_list

fastq_path <- fastqs[1]
generate_subsetted_fastqs(fastq_path, all_classed, unclassified_reads_out,
                          outdir_organism=organism_dir,
                          outdir_kingdom=kingdom_dir )
