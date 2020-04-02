#
#
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(data.table)} )
  
kraken_files <- snakemake@input[['krakenfile']]
output_file <-  snakemake@output[['kraken_classified_filtered']]
classified_count_file <-  snakemake@output[['read_count']]
kmer_len <- snakemake@params[['kmer_len']]

read_and_filter_kraken_file <- function(file_name, kmer_len) {

  df <-  data.table::fread(file = file_name) %>% 
    as_tibble() %>% 
    setNames(c("classified","seq_id","tax_id","seq_length","score", "kmer_LCA") ) 
  
  df %>% 
    filter(classified=="C") %>% 
    mutate(score=as.double(str_remove(score,"P=")),
      matches=(seq_length-kmer_len)*score+0.5, #convert score to number of matches.
      type=str_remove(basename(file_name),"\\.txt$") ) %>% 
    select(seq_id, tax_id, matches, classified, type)
}

kraken_files %>%
  map_dfr(~{read_and_filter_kraken_file(.x,kmer_len)}) %T>%  
  print() %>% 
  write_tsv(path = output_file) %>% 
  group_by(type) %>% 
  summarise(count=n()) %>% 
  print() %>%
  write_tsv(path = classified_count_file)



# kraken_files <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/kraken_filtered_viruses.txt"
# output_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/kraken_filtered_viruses2.txt"
# classified_count_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/count.txt"
# kmer_len <- 30
# 
# kraken_es <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/kraken_filtered_cows.txt"
# kraken_files <- c(kraken_files, kraken_es)