# Author: Pernilla Ericsson (pernilla.ericsson@gu.se, clinicalgenomics@gu.se)
# Date: 2020-12-22

suppressPackageStartupMessages({
  library(tidyverse);library(magrittr);library(data.table);library(ShortRead)
} )

read_count <- snakemake@input[['read_count']]
trimmed_read_count <- snakemake@input[['trimmed_read_count']]

mincount <- snakemake@params[['mincount']]

tableview_out <- snakemake@output[['tableview']]
classified_reads_mincount <- snakemake@output[['classified_reads_mincount']]
read_count_out <- snakemake@output[['read_count']]

organism_dir <- snakemake@output[['organism_dir']]
kingdom_dir <- snakemake@output[['kingdom_dir']]


quit_when_empty <- function(tableview, classified_reads,
                            read_count, organism_outdir, kingdom_outdir){
  # Function to handle if there are no classified reads.
  if(!dir.exists(organism_outdir)) dir.create(organism_outdir, recursive = TRUE)
  if(!dir.exists(kingdom_outdir)) dir.create(kingdom_outdir, recursive = TRUE)
  write_tsv(tibble(), tableview)
  write_tsv(tibble(), classified_reads)
  write_tsv(tibble(count=0), read_count)
  quit(save = "no", status = 0)
}

# Check if read count is empty and quit.
if (file.size(read_count) == 0) {
  quit_when_empty(tableview=tableview_out,
                  classified_reads=classified_reads_mincount,
                  read_count=read_count_out, 
                  organism_outdir=organism_dir, 
                  kingdom_outdir=kingdom_dir)
}

read_total <- read_tsv(trimmed_read_count) %>% pull(count)

# Read in the classified reads dataframe with readcount included.
read_count_df <-data.table::fread(file = read_count, fill=TRUE, sep = "\t",
                                  header=TRUE) %>% as_tibble()

# Add a column with percent of total reads.
tableview <- 
  read_count_df %>% 
  mutate(superkingdom=str_trim(superkingdom) %>% str_replace_all(" ", "_") ) %>% 
  mutate(superkingdom=tolower(superkingdom)) %>%
  mutate(tax_id=as.character(tax_id)) %>% 
  mutate(across(where(is.character), ~replace_na(., "NA") )) %>% 
  group_by(superkingdom, organism, tax_id) %>% 
  summarize(read_count=sum(read_count, na.rm = TRUE), .groups = "drop")  %>% 
  mutate(percent=round((read_count/read_total)*100, digits=3), .before=read_count) %>% 
  filter(read_count>mincount) 

if(nrow(tableview)==0){
  quit_when_empty(tableview=tableview_out,
                  classified_reads=classified_reads_mincount,
                  read_count=read_count_out, 
                  organism_outdir=organism_dir, 
                  kingdom_outdir=kingdom_dir)
} else {
  if(!dir.exists(organism_dir)) dir.create(organism_dir, recursive = TRUE)
  if(!dir.exists(kingdom_dir)) dir.create(kingdom_dir, recursive = TRUE)
}

tableview %>% 
  write_tsv(tableview_out, col_names = TRUE)

tableview %>% select(read_count) %>%
  summarize(count=sum(read_count)) %>% 
  write_tsv(read_count_out, col_names = TRUE)

# Find read ids at a taxid with number of reads higher than mincount and write to a file.
unique_tax_id <- tableview %>% pull(tax_id)
classified_df <-
  read_count_df %>% select(tax_id, seq_id, superkingdom, organism) %>% 
  filter(tax_id %in% unique_tax_id) %>% 
  write_tsv(classified_reads_mincount, col_names = TRUE)

# Create filed called after taxid containing read names
organism_df_list <- 
  classified_df %$% 
  split(., tax_id)
imap(organism_df_list, ~{
    outfile_organism <- glue::glue("{organism_dir}/organism_{.y}.tsv")
    write_tsv(.x, path = outfile_organism, col_names = TRUE)
  })

# Create filed called after kingdom containing read names
kingdom_df_list <- 
  classified_df %$% 
  split(., superkingdom)
imap(kingdom_df_list, ~{
  outfile_kingdom <- glue::glue("{kingdom_dir}/kingdom_{.y}.tsv")
  write_tsv(.x, path = outfile_kingdom, col_names = TRUE)
})


# trimmed_read_count <-"/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1/snakemake_results_sample_1/stats_PE_RNA/stage1/trimming/count_bbduk_unmerged_reads_trimmed_raw.txt"
# 
# read_count <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1/snakemake_results_sample_1/PE_RNA/stage8/readcount.tsv"

# mincount <- 4
# 
# tableview_out <- "/Users/pernillaericsson/Desktop/parca_tests/tableview.tsv"
# classified_reads_mincount <- "/Users/pernillaericsson/Desktop/parca_tests/classified_reads_mincount.tsv"
# read_count_out <- "/Users/pernillaericsson/Desktop/parca_tests/classified_reads_count.tsv"
# organism_dir <- "/Users/pernillaericsson/Desktop/parca_tests/organism_dir"
# kingdom_dir <- "/Users/pernillaericsson/Desktop/parca_tests/kingdom_dir"



# trimmed_read_count_unmerged <-"/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1/snakemake_results_sample_1/stats_PE_RNA/stage1/trimming/count_bbduk_unmerged_reads_trimmed_raw.txt"
# trimmed_read_count_merged <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1/snakemake_results_sample_1/stats_PE_RNA/stage1/trimming/count_bbduk_merged_reads_trimmed.txt"

# Read in the read count stats with regard to PE or SE
# if (SE_or_PE == "SE") {
#   trimmed_read_count <- snakemake@input[['trimmed_read_count']]
#   
#   read_total <- read_tsv(trimmed_read_count) %>% pull(count)
#   
# } else if (SE_or_PE == "PE"){
#   trimmed_read_count_unmerged <- snakemake@input[['trimmed_read_count_unmerged']]
#   trimmed_read_count_merged <- snakemake@input[['trimmed_read_count_merged']]
#   
#   read_total <- 
#     c(trimmed_read_count_unmerged, trimmed_read_count_merged) %>% 
#     map_dfr(~read_tsv(.x)) %>% select(count) %>% colSums() %>% unname()
#   
# } else {
#   quit(save="no", status=0)
# }
