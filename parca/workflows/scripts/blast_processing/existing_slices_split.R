

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
} )

# higher_file <-"/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage4/genusspeciessplit/above_species_classed.txt"
# min_tax_id_count <- 2
# existing_slice_list <- read_tsv("/Users/pernillaericsson/Documents/medair1/home/xerpey/gislices", col_names = "files") %>% pull(files)
# missing_slice_file <- "/Users/pernillaericsson/missing.txt"
# detected_slice_file<- "/Users/pernillaericsson/detected.txt"

higher_file <- snakemake@input[["higher"]]

existing_slice_path <- snakemake@params[["existing_slice_path"]]
min_tax_id_count <- snakemake@params[["min_tax_id_count"]]

detected_slice_file <- snakemake@output[["detected"]]
missing_slice_file <- snakemake@output[["missing"]]
count_file <- snakemake@output[["count"]]

existing_slice_list <- list.files(path=existing_slice_path)
df_rank_kingdom_genus_species <- read_tsv(higher_file) 

detected <- 
  df_rank_kingdom_genus_species %>% 
  filter(tax_id %in% existing_slice_list) %>% 
  pull(tax_id) %>% unique() %>% 
  tibble::enframe(name = NULL) %>% 
  write_tsv(detected_slice_file,col_names = FALSE)

missing <- 
  df_rank_kingdom_genus_species %>%
  filter(!tax_id %in% existing_slice_list) %>% 
  group_by(tax_id) %>% 
  summarise(count=n()) %>% 
  filter(count>min_tax_id_count) %>%
  pull(tax_id) %>% unique() %>% tibble::enframe(name = NULL) %>%  
  write_tsv(missing_slice_file,col_names = FALSE)


df_count <- 
  tibble(type=c("detected","missing"),count=c(nrow(detected),nrow(missing)) ) %>% 
  bind_rows(tibble(type="total",count=sum(c(nrow(detected),nrow(missing))))) %>% 
  write_tsv(count_file)
