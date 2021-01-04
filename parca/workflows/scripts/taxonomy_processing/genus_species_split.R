#
#
# Author: Pernilla Ericsson (pernilla.ericsson@gu.se)
# Date: 2020-05

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(data.table) 
} )

# combined_doublets_singletons_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage4/taxonomy_processing/combined_doublets_singletons.txt"
# 
# species_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage4/genusspeciessplit/species_classed.txt"
# higher_file <-"/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage4/genusspeciessplit/above_species_classed.txt"
# read_count_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage4/genusspeciessplit/read_count_species_genus_higher.txt"

# combined_doublets_singletons_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/200819_demo/snakemake_results_SRR1761912/PE_RNA/stage4/taxonomy_processing/combined_doublets_singletons.txt"
# species_file <- "/Users/pernillaericsson/Desktop/species_classed.txt"
# higher_file < -"/Users/pernillaericsson/Desktop/above_species_classed.txt"
# read_count_file <- "/Users/pernillaericsson/Desktop/read_count_species_genus_higher.txt"

combined_doublets_singletons_file <- snakemake@input[["combined_doublets_singletons"]]

species_file <- snakemake@output[["species"]]
higher_file <- snakemake@output[["higher"]]
read_count_file <- snakemake@output[["read_count"]]

df <- data.table::fread(file = combined_doublets_singletons_file, 
                        fill=TRUE, 
                        sep = "\t",
                        header=TRUE) %>% 
  as_tibble() 

if(nrow(df)==0) {
  write_tsv(tibble(), species_file)
  write_tsv(tibble(), higher_file)
  write_tsv(tibble(type="Total", count=0), read_count_file)
  quit(save = "no", status = 0)
}

df %<>% 
  filter(
  tax_id != 1,
  classified == "C" ) %>% 
  mutate(
    species = ifelse(order %in% "Primates", "Homo_sapiens", species),
    tax_id = ifelse(order %in% "Primates", 9606, tax_id) )

df %<>% 
  mutate_if(is.character, str_trim) %>% 
  mutate(
    species=str_replace_all(species, " ", "_"),
    genus=str_replace_all( genus, " ", "_"),
    family=str_replace_all( family, " ", "_"),
    order=str_replace_all( order, " ", "_"),
    class=str_replace_all( class, " ", "_"),
    superkingdom=str_replace_all( superkingdom, " ", "_"),
    ) 

df_tmp_filter <- df %>% select(species, genus, family, order, class, superkingdom)

df <- df[rowSums(is.na(df_tmp_filter)) != ncol(df_tmp_filter), ]

df

df_rank <- 
  df %>% 
  rowwise() %>%
  mutate(rank = pmap_chr(.l = list(species, genus, family, order, class, superkingdom),
                         .f = ~ first(names(discard(c(s=..1, g=..2 , f=..3, o=..4 , c=..5, k=..6), ~is.na(.x) ) ) ) ),
         organism = pmap_chr(.l = list(species, genus, family, order, class, superkingdom),
                             .f = ~ first(discard(c(s=..1, g=..2 , f=..3, o=..4 , c=..5, k=..6), ~is.na(.x) ) ) )  ) %>% 
  ungroup()

df_rank_kingdom_genus_species <- 
  df_rank %>% 
  mutate(
    higher_genus_species = ifelse(
      rank %in% "s", "species", ifelse(rank %in% "g", "genus", "higher" ) ),
    higher_species = ifelse(
      rank %in% c("g", "f", "o", "c", "k"), "higher", "species"   ) )

total <- 
  df_rank_kingdom_genus_species %>% 
  nrow() %>%
  bind_cols(type="total",count=.)

count_higher_genus_species <- 
  df_rank_kingdom_genus_species %>% 
  group_by(higher_genus_species) %>% 
  summarise(count=n()) %>% 
  rename("type"=higher_genus_species) %>% 
  arrange(match(type, c("species", "genus", "higher")) ) %>% 
  bind_rows(total) %>% 
  write_tsv(path=read_count_file)


df_rank_kingdom_genus_species %>% 
  filter( higher_species %in% "species") %>% 
  mutate(
    score=1,
    type="kmer") %>% 
  select(seq_id, tax_id, score, rank, organism, type) %>% 
  write_tsv(path = species_file)

df_rank_kingdom_genus_species %>% 
  filter( higher_species %in% "higher") %>% 
  mutate(
    score=1,
    type="kmer") %>% 
  select(seq_id, tax_id, score, rank, organism, type) %>% 
  write_tsv(path = higher_file)

#df %>% head() %>% dput()

