#
#

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(data.table) } )

kraken_file <- snakemake@input[['kraken']]
kaiju_file <- snakemake@input[['kaiju']]

kraken_doublets_file <- snakemake@output[['kraken_doublets']]
kaiju_doublets_file <- snakemake@output[['kaiju_doublets']]
singletons_file <- snakemake@output[['singletons']]
merged_file <- snakemake@output[['merged']]
doublet_count_file <- snakemake@output[['read_count_doublet']]
singletons_count_file <- snakemake@output[['read_count_singletons']]


kmer_files <- c(kraken_file, kaiju_file)
#print(kmer_files)

read_files <- function(file_name) {
  #print(file_name)
  df <- data.table::fread(file = file_name, 
                               fill=TRUE, 
                               sep = "\t" ) %>% 
    as_tibble() %>% 
    setNames(c("classified", "seq_id", "tax_id", "matches", "kmer_counter", "type") )
  }


merged_df <- kmer_files %>%
  map_dfr(~{read_files(.x)}) %>% 
  write_tsv(path=merged_file, col_names = FALSE) 
  #print()

doublet <-
  merged_df %>% 
  group_by(seq_id) %>% 
  summarise(hits=n()) %>% 
  filter(hits==2) %>%
  pull(seq_id)
write_tsv(tibble("Doublets", length(doublet)), path = doublet_count_file, col_names = FALSE)


merged_df_doublets <- 
  merged_df %>% 
  filter(seq_id %in% doublet) %>% 
  arrange(seq_id)

kraken_doublets <- 
  merged_df_doublets %>% 
  filter(kmer_counter=="kraken") %>% 
  write_tsv(path=kraken_doublets_file, col_names = FALSE)
  
kaiju_doublets <- 
  merged_df_doublets %>% 
  filter(kmer_counter=="kaiju") %>% 
  write_tsv(path=kaiju_doublets_file, col_names = FALSE)

merged_df_singletons <- 
  merged_df %>% 
  filter(!(seq_id %in% doublet)) %>% 
  arrange(seq_id) %>% 
  write_tsv(path=singletons_file, col_names = FALSE) %>% 
  group_by(kmer_counter) %>% 
  summarise(count=n()) 

merged_df_singletons %>% 
  summarize("count"=sum(count))  %>% 
  mutate(kmer_counter="Total") %>% 
  bind_rows(merged_df_singletons,.) %>% 
  #print() %>% 
  write_tsv(path = singletons_count_file)


# kraken_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage3/kraken/kraken_filtered_classified.txt"
# kaiju_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage3/kaiju/kaiju_filtered_classified.txt"
# kraken_doublets_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage3/kraken_dbl_test.txt"
# kaiju_doublets_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage3/kaiju_dbl_test.txt"
# singletons_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage3/singletons_test.txt"
# merged_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage3/merged_test.txt"
# doublet_count_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage3/count_doublets_test.txt"
# singletons_count_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage3/count_singletons_test.txt"