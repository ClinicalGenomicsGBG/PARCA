# Author: Pernilla Ericsson (pernilla.ericsson@gu.se, clinicalgenomics@gu.se)
# Date: 2020-12-18

suppressPackageStartupMessages({
  library(tidyverse);library(magrittr);library(data.table)
} )

# read_count <- "/Users/pernillaericsson/Desktop/readcount.tsv"
# read_total <- 1104870
# mincount <- 4
# organism_dfs_dir <- "/Users/pernillaericsson/Desktop/test_df_out"


read_count <- snakemake@input[['read_count']]

SE_or_PE <- snakemake@params[['SE_or_PE']]
mincount <- snakemake@params[['mincount']]

tableview_out <- snakemake@output[['tableview']]
organism_dfs_dir <- snakemake@output[['organism_dfs_dir']]


if (SE_or_PE == "SE") {
  trimmed_read_count <- snakemake@input[['trimmed_read_count']]
  read_total <- read_tsv(trimmed_read_count) %>% pull(count)

} else if (SE_or_PE == "PE"){
  trimmed_read_count_unmerged <- snakemake@input[['trimmed_read_count_unmerged']]
  trimmed_read_count_merged <- snakemake@input[['trimmed_read_count_merged']]

  read_total <- 
    c(trimmed_read_count_unmerged, trimmed_read_count_merged) %>% 
    map_dfr(~read_tsv(.x)) %>% select(count) %>% colSums() %>% unname()

} else {
  quit(save="no", status=0)
}

read_count_df <-data.table::fread(file = read_count, fill=TRUE, sep = "\t",
                                  header=TRUE) %>% as_tibble()

if(nrow(read_count_df)==0){
  if(!dir.exists(organism_dfs_dir)) dir.create(organism_dfs_dir, recursive = TRUE)
  write_tsv(tibble(), tableview_out)
  quit(save = "no", status = 0)
}

tableview <- 
  read_count_df %>% 
  group_by(superkingdom, organism, tax_id) %>% 
  summarize(read_count=sum(read_count) )  %>% 
  mutate(percent=(read_count/read_total)*100) %>% 
  relocate(read_count, .after = last_col()) %>% 
  filter(read_count>mincount)

unique_tax_id <- tableview %>% pull(tax_id)

organism_df_list <-
  read_count_df %>% select(tax_id, seq_id, superkingdom, organism) %>% 
  filter(tax_id %in% unique_tax_id) %$% split(., tax_id )

tableview %>% 
  write_tsv(tableview_out, col_names = TRUE)

# Function to save the df with the name of the tax_id
save_by_taxid <- function(organism_df, outfile){
  out_name <- paste0(organism_dfs_dir,"/",outfile)
  write_tsv(organism_df, out_name, col_names = TRUE) 
}

if(!dir.exists(organism_dfs_dir)) dir.create(organism_dfs_dir, recursive = TRUE)

map2(organism_df_list, names(organism_df_list), ~save_by_taxid(.x, .y) )
