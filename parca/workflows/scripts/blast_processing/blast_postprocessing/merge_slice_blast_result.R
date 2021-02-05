# Script for readling multiple blast result files and grouping the reads and taxids in order to sum their bitscore and returning the read and taxid with the best score
# Author: Pernilla Ericsson
# Date: 2020-05

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(data.table) 
} )

files <- snakemake@input[['blast_output']]
best_blast <- snakemake@output[['best_blast']]


# Function to read a blast result file
read_blast_results <- function(file_name) {
  df <- data.table::fread(file_name)
  if (nrow(df)==0){
    df <- tibble()
  } else {
  df %<>% 
    as_tibble() %>% 
    setNames(c("qseqid", "stitle", "evalue", "length", "nident", "sseqid", "bitscore", "staxids") ) %>% 
    mutate(staxids=str_remove(staxids, ".*;")) 
  }
}

# Read all blast result files and remove duplicate rows using distinct().
df_all <- 
  files %>% 
  map_dfr(~read_blast_results(.x)) %>% 
  distinct() 

# If no blast results are found, write an empty table to the output name and quit the program.
if (nrow(df_all)==0) {
  write_tsv(tibble(),best_blast)
  quit(save = "no", status = 0)
}

# Sum all bitscores for reads assigned to the same taxid multiple times.
df_bitscore <- 
  df_all %>%
  select(qseqid, staxids, bitscore) %>% 
  group_by(qseqid, staxids) %>% 
  summarize(score=sum(bitscore)) %>%
  ungroup() 

# select the read for which the taxid has the highest score.
df_bitscore_best <- 
  df_bitscore %>% 
  group_by(qseqid) %>% 
  arrange(desc(score)) %>%
  dplyr::slice(1)

write_tsv(df_bitscore_best,best_blast)


# file_name = "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage6/sliceblastout/9984__1"
# file_name0 = "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage6/sliceblastout/9984__0"
# file_name1="/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage6/sliceblastout/9903__0"
# 
# files <- c(file_name, file_name0,file_name1)
#best_blast <- "/Users/pernillaericsson/Desktop/a.tsv"