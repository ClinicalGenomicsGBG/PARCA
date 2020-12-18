# Author: Pernilla Ericsson (pernilla.ericsson@gu.se, clinicalgenomics@gu.se)
# Date: 2020-12-18

suppressPackageStartupMessages({
  library(tidyverse);library(magrittr);library(data.table)
} )

# read_count <- "/Users/pernillaericsson/Desktop/readcount.tsv"
# read_total <- 1104870

read_count <- snakemake@input[['read_count_df']]

SE_or_PE <- snakemake@params[['SE_or_PE']]
mincount <- snakemake@params[['mincount']]

tableview_out <- snakemake@output[['tableview']]
organism_dfs_dir <- snakemake@output[['organism_dfs_dir']]

if (SE_or_PE == "SE") {
  trimmed_read_count <- snakemake@input[['trimmed_read_count']]
  read_total <- read_tsv(trimmed_read_count) %>% pull(count)
  
} elif (SE_or_PE == "PE"){
  trimmed_read_count_unmerged <- snakemake@input[['trimmed_read_count_unmerged']]
  trimmed_read_count_merged <- snakemake@input[['trimmed_read_count_merged']]
  
  read_total <- 
    c(trimmed_read_count_unmerged, trimmed_read_count_merged) %>% 
    map_dfr(~read_tsv(.x)) %>% select(count) %>% colSums() %>% unname()
} else {
  quit(save="no", status=0)
}

read_count_df <-data.table::fread(file = read_count, 
                  fill=TRUE, 
                  sep = "\t",
                  header=TRUE) %>% as_tibble()

tableview <- 
  read_count_df %>% 
  group_by(superkingdom, organism, tax_id, tax_id_sgft_lower) %>% 
  summarize(read_count=sum(read_count) )  %>% 
  mutate(percent=(read_count/read_total)*100) %>% 
  relocate(read_count, .after = last_col()) %>% 
  filter(read_count>mincount)

unique_tax_id <- tableview %>% pull(tax_id_sgft_lower)

organism_df_list <-
  read_count_df %>% select(seq_id, superkingdom ,tax_id_sgft_lower) %>% 
  filter(tax_id_sgft_lower %in% unique_tax_id) %$% split(., tax_id_sgft_lower )

tableview %>% 
  write_tsv(tableview_out, col_names = TRUE)