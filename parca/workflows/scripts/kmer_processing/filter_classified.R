#
#
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(data.table)} )
  
kmer_files <- snakemake@input[['files']]
output_file <-  snakemake@output[['classified_filtered']]
classified_count_file <-  snakemake@output[['read_count']]
selected_program <- snakemake@params[['program']]
selected_kmer_len <- 30

read_and_filter_file <- function(file_name, kmer_len=30, program="") {
  
  if(program=="kraken"){
  df <-  data.table::fread(file = file_name, sep = "\t") %>% 
    as_tibble() %>% 
    setNames(c("classified","seq_id","tax_id","seq_length","score", "kmer_LCA") ) 
  df %>% 
    filter(classified=="C") %>% 
    mutate(score=as.double(str_remove(score,"P=")),
      matches=(seq_length-kmer_len)*score+0.5, #convert score to number of matches.
      kmer_counter="kraken",
      type=str_remove(basename(file_name),"\\.txt$") ) %>% 
    select(seq_id, tax_id, matches, classified,kmer_counter, type) } else if (program=="kaiju") {
      
    df <- data.table::fread(file = kaiju_file, fill=TRUE, sep = "\t") %>% 
      as_tibble() %>% 
      setNames(c("classified","seq_id","tax_id","seq_length_score","tax_ids_best","accession_best","fragment_seq") ) 
    df %>% 
      filter(classified=="C") %>% 
      separate(fragment_seq,into=c("best_seq","other_seqs"),extra = "merge", fill="right", sep = ",") %>% 
      mutate(
        prot_length=nchar(best_seq),
        matches=nchar(best_seq)*3,
        kmer_counter="kaiju",
        type=str_remove(basename(file_name),"\\.txt$") ) %>% 
      select(seq_id, tax_id, matches, classified,kmer_counter, type)
    }
  
}

best_classifications <- 
  kmer_files %>%
  map_dfr(~{read_and_filter_file(.x,kmer_len=selected_kmer_len,program=selected_program)}) %>%
  group_by(seq_id) %>% arrange(desc(matches)) %>% slice(1) %>% ungroup() %>% 
  write_tsv(path = output_file) %>% 
  print() %>%
  group_by(type) %>%
  summarise(count=n()) 

best_classifications %>% 
  summarize("count"=sum(count))  %>% 
  mutate(type="Total") %>% 
  bind_rows(best_classifications,.) %>% 
  print() %>% 
  write_tsv(path = classified_count_file)



# kaiju_file="/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage3/kaiju/kaijuresults_progenomes_complete_80_18.txt"
# selected_program="kraken"
# 
# kraken_files <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/kraken_filtered_viruses.txt"
# output_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/kraken_filtered_viruses2.txt"
# classified_count_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/count.txt"
# kmer_files <- kraken_files

# kraken_es <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/kraken_filtered_cows.txt"
# kraken_files <- c(kraken_files, kraken_es)

