#
# Script for separating unclassed species, genus or family, i.e SGF, from combined kaiju and kraken result. The reads with unclassed SGF are placed in the singletons file.
# Author: Pernilla Ericsson (pernilla.ericsson@gu.se)
# Date: 2020-05

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(data.table) 
  } )

unfiltered_file <- snakemake@input[['combined_unfiltered']]
kraken_doublets_file <- snakemake@input[['kraken_doublets']]
kaiju_doublets_file <- snakemake@input[['kaiju_doublets']]
singletons_file <- snakemake@input[['singletons']]

filtered_SGF_empty_removed_file <- snakemake@output[['combined_SGF_empty_filter']]
singletons_SGF_empty_file <- snakemake@output[['singletons_added_SGF_empty']]

## Combined kaiju and kraken
df <- data.table::fread(file = unfiltered_file, 
                        fill=TRUE, 
                        sep = "\t" ) %>% 
  as_tibble() %>% 
  setNames(c("classified", "seq_id", "tax_id", "tax_names_lineage") )

if (nrow(df)==0) {
  write_tsv(tibble(),filtered_SGF_empty_removed_file)
  df_SGF_empty_only <- tibble()
  
} else {
  
  df_clean <- 
    df %>% 
    separate(col=tax_names_lineage, 
                  into=c("superkingdom","class","order","family","genus","species","other"),
                  sep=";",
                  extra = "merge") %>% select(-other) %>% 
    mutate_if(is.character, str_trim) %>% mutate_all(na_if,"NA")
  
  df_SGF_empty_only <- 
    df_clean %>% 
    filter(
      is.na(species),
      is.na(genus),
      is.na(family) )
  
  df_SGF_empty_removed <- 
    df_clean %>% 
    filter(
      !is.na(species),
      !is.na(genus),
      !is.na(family)
    ) %>% write_tsv(filtered_SGF_empty_removed_file)
}

## Add reads with no names for species, genus and family to singletons file
read_files <- function(file_name) {
  #print(file_name)
  df <- data.table::fread(file = file_name, 
                          fill=TRUE, 
                          sep = "\t" ) %>% 
    as_tibble() %>% 
    setNames(c("classified", "seq_id", "tax_id", "matches", "kmer_counter", "type") )
}

df_doublets <- c(kraken_doublets_file,kaiju_doublets_file) %>%
  map_dfr(~{read_files(.x)})

df_singletons <- 
  read_files(singletons_file) 


if (nrow(df_singletons)==0) {
  df_singletons <- tibble()
} else{
  df_singletons %<>% 
    mutate(tax_id=as.character(tax_id) )
}


if (nrow(df_SGF_empty_only)==0) {
  vec_SGF_empty_only <- c()
} else {
  vec_SGF_empty_only <- df_SGF_empty_only %>% pull(seq_id)
}

if (nrow(df_doublets)==0) {
  write_tsv(tibble(),singletons_SGF_empty_file)
} else {
  df_SGF_empty_best_singletons <-
    df_doublets %>%
    filter(seq_id %in% vec_SGF_empty_only ) %>%
    arrange(seq_id, desc(matches)) %>%
    group_by(seq_id) %>%
    slice(1) %>% ungroup() %>% 
    mutate(tax_id=as.character(tax_id))
  
  df_singletons_SGF_empty <- 
    bind_rows(df_SGF_empty_best_singletons, df_singletons) %>% 
    write_tsv(singletons_SGF_empty_file)
}


# unfiltered_file <-"/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage4/comparison/combined_kraken_kaiju_names_unfiltered.txt"
# singletons_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage4/comparison/singletons.txt"
# kraken_doublets_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage4/comparison/kraken_doublets.txt"
# kaiju_doublets_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage4/comparison/kaiju_doublets.txt"
# 
# 
# filtered_SGF_empty_removed_file <-"/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage4/filter_combined/combined_kraken_kaiju_names.txt"
# singletons_SGF_empty_file <-"/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage4/filter_combined/combined_kraken_kaiju_names.txt"

# vec_singletons_SGF_empty_taxid <-
#   df_singletons_SGF_empty %>%
#   pull(tax_id) %>% unique()

# 
# df_singletons_SGF_empty_names <-
#   getTaxonomy(vec_singletons_SGF_empty_taxid,
#                  sql_db_file,
#                  desiredTaxa = c("superkingdom", "phylum", "class", "order", "family",
#                                  "genus", "species") ) #superkingdom,class,order,family,genus,species
# 
# SGF_empty_singletons_names_tax_ids <-
#   rownames(df_singletons_SGF_empty_names) %>% str_trim()
# 
# df_tax_names <-
#   df_singletons_SGF_empty_names %>%
#   as_tibble() %>%
#   bind_cols(tax_id=SGF_empty_singletons_names_tax_ids,.)
# 
# # df_tax_names %>% filter(is.na(superkingdom) )
# singletons_SGF_empty <- 
#   left_join(df_singletons_SGF_empty,df_tax_names, by="tax_id") %>% 
#   write_tsv(singletons_SGF_empty_file)
# 
# # singletons_SGF_empty %>% 
# #   filter(seq_id == "k141_11")
# 
# singletons_SGF_empty %>% 
#   unite("lowest_vec", c(tax_id, superkingdom, phylum, class, order, family, genus, species), na.rm = TRUE, remove = FALSE, sep="__") %>% 
#   mutate(lowest = str_remove(lowest_vec,"\\w+__") ) %>% 
#   mutate(final_tax_id=lowest %>% getId(., sqlFile = sql_db_file, onlyScientific = TRUE) ) #%>% 
#   #filter(str_detect(final_tax_id, ",") )
#          