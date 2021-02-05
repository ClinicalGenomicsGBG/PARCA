#
#
# Author: Pernilla Ericsson (pernilla.ericsson@gu.se)
# Date: 2020-05-13

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(data.table) 
  library(dplyr)
} )

all_classed <- snakemake@input[['all_classed']]

all_classed_read_taxid <- snakemake@output[['all_classed_read_taxid']]
type_summary <- snakemake@output[['type_summary']]

# Read the file.
df_all <- data.table::fread(file = all_classed, 
                              fill=TRUE, 
                              sep = "\t",
                              header=TRUE) %>% 
  as_tibble()


df_all %>% 
  filter(!tax_id %in% 0) %>% 
  mutate(classified="C",
         taxid_score=paste(tax_id, score, sep = ":") ) %>% 
  select(classified, seq_id, tax_id, score, taxid_score) %>% 
  write_tsv(all_classed_read_taxid, col_names = FALSE)

df_all %>% 
  filter(!tax_id %in% 0) %>% 
  group_by(type) %>% tally() %>% 
  bind_rows(summarise(., 
                      across(where(is.numeric), sum),
                      across(where(is.character), ~"Total"))) %>% 
  write_tsv(type_summary, col_names = TRUE) 

# all_classed <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage7/kmer_species_subsetblast_blastnt_classed.txt"
# all_classed_read_taxid <- "/Users/pernillaericsson/Desktop/all_classed.txt"