#
#

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(data.table) 
} )

singletons_genus_names <- snakemake@input[['singletons_genus_names']]

singletons_genus_names_reformat <- snakemake@output[['singletons_genus_names_reformat']]

#singletons_genus_names <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage4/taxonomy_processing/singletons_genus_names.txt"

df <- data.table::fread(file = singletons_genus_names, 
                        fill=TRUE, 
                        sep = "\t" ) %>% 
  as_tibble() %>% 
  setNames(c("classified", "seq_id", "tax_id", "tax_names_lineage") )


df_clean <- 
  df %>% 
  separate(col=tax_names_lineage, 
           into=c("superkingdom","class","order","family","genus","species","other"),
           sep=";",
           extra = "merge") %>% select(-other) %>% 
  mutate_if(is.character, str_trim) %>% mutate_all(na_if,"NA") %>% 
  write_tsv(singletons_genus_names_reformat)

