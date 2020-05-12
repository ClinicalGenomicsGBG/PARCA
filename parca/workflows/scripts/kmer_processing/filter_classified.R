# Script for processing kraken or kaiju result and selecting the best classification of a read by adding a score called matches and selecting the highest one.
# Author: Pernilla Ericsson (pernilla.ericsson@gu.se)
# Date: 2020-05

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr) } )
  #library(data.table)} )
  
kmer_files <- snakemake@input[['files']]
output_file <-  snakemake@output[['classified_filtered']]
classified_count_file <-  snakemake@output[['read_count']]
selected_program <- snakemake@params[['program']]
selected_kmer_len <- 30


read_and_filter_file <- function(file_name, kmer_len=30, program="") {
  #print(file_name)
  if(program=="kraken"){
    df <-  read_tsv(file = file_name, 
                             col_names = c("classified",
                                           "seq_id",
                                           "tax_id",
                                           "seq_length",
                                           "score", 
                                           "kmer_LCA")) %>% 
      as_tibble() 
    
    df %<>% 
      filter(classified=="C")
    
    if (nrow(df)==0) {
      return(df)
    }
    
    df %>% 
      mutate(score=as.double(str_remove(score,"P=")),
        matches=(seq_length-kmer_len)*score+0.5, #convert score to number of matches.
        kmer_counter="kraken",
        type=str_remove(basename(file_name),"\\.txt$") ) %>% 
      select(classified, seq_id, tax_id, matches, kmer_counter, type) 
    
  } else if (program=="kaiju") {
      
      df <- read_tsv(
        file = file_name, 
        col_names=c("classified",
                    "seq_id",
                    "tax_id",
                    "seq_length_score",
                    "tax_ids_best",
                    "accession_best",
                    "fragment_seq") 
      ) %>% 
        as_tibble()
      
      df %<>% 
        filter(classified=="C") 
      
      if (nrow(df)==0) {
        return(df)
      }
      
      df %>% 
        separate(fragment_seq,into=c("best_seq","other_seqs"),extra = "merge", fill="right", sep = ",") %>% 
        mutate(
          prot_length=nchar(best_seq),
          matches=nchar(best_seq)*3,
          kmer_counter="kaiju",
          type=str_remove(basename(file_name),"\\.txt$") ) %>% 
        select(classified, seq_id, tax_id, matches, kmer_counter, type)
    }
  
}

best_classifications <- 
  kmer_files %>%
  map_dfr(~{read_and_filter_file(.x,kmer_len=selected_kmer_len,program=selected_program)}) 

if(nrow(best_classifications)==0) {
  write_tsv(best_classifications,path = output_file)
  stats_df <- tibble(type="Total", count=0)
  write_tsv(stats_df,path = classified_count_file)
  quit(save="no", status=0)
}

best_classifications %<>%
  group_by(seq_id) %>% arrange(desc(matches)) %>% slice(1) %>% ungroup() %>% 
  write_tsv(path = output_file) %>% 
  #print() %>%
  group_by(type) %>%
  summarise(count=n()) 

best_classifications %>% 
  summarize("count"=sum(count))  %>% 
  mutate(type="Total") %>% 
  bind_rows(best_classifications,.) %>% 
  #print() %>% 
  write_tsv(path = classified_count_file)



# 
# file_name <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage3/kaiju/kaijuresults_genbank_viruses_171017_kaiju_db_75_15.txt"
# 
# df <- read_tsv(
#   file = file_name, 
#   col_names=c("classified",
#               "seq_id",
#               "tax_id",
#               "seq_length_score",
#               "tax_ids_best",
#               "accession_best",
#               "fragment_seq") 
# ) %>% 
#   as_tibble()


# kaiju_file="/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage3/kaiju/kaijuresults_progenomes_complete_80_18.txt"
# selected_program="kraken"
#  
#  
#file_name <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage3/kaiju/kaijuresults_genbank_eukaryotes_171014_kaiju_db_85_20.txt"
# kraken_files <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/kraken_filtered_viruses.txt"
# output_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/kraken_filtered_viruses2.txt"
# classified_count_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/count.txt"
# kmer_files <- kraken_files

# kraken_es <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/kraken_filtered_cows.txt"
# kraken_files <- c(kraken_files, kraken_es)
# 
#     # df <- data.table::fread(file = file_name, 
#                         fill=TRUE, 
#                         sep = "\t",
#                         col.names=c("classified",
#                                    "seq_id",
#                                    "tax_id",
#                                    "seq_length_score",
#                                    "tax_ids_best",
#                                    "accession_best",
#                                    "fragment_seq") ) %>% 
# as_tibble() %>%
# setNames(c("classified","seq_id","tax_id","seq_length_score","tax_ids_best","accession_best","fragment_seq") )
# 
# 
# kmer_files <- "/Users/pernillaericsson/Desktop/empty_test.txt"
# output_file <-  "/Users/pernillaericsson/Desktop/empty_test_filtered.txt"
# classified_count_file <-  "/Users/pernillaericsson/Desktop/count_empty_test.txt"
# selected_program <- "kraken"


