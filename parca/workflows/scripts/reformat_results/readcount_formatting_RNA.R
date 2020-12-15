# Author: Pernilla Ericsson (pernilla.ericsson@gu.se, clinicalgenomics@gu.se)
# Date: 2020-12-14

suppressPackageStartupMessages({
  library(tidyverse);library(magrittr);library(data.table)
} )

classed_reads_path <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/200819_demo/snakemake_results_a/SE_RNA/stage8/all_classed_read_taxid_names.txt"

coverage_path <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/200819_demo/snakemake_results_a/SE_RNA/stage2/pileup/bbmap_cov.txt"

outfile <- "/Users/pernillaericsson/Desktop/readcount.tsv"

df_classed_reads <- data.table::fread(file = classed_reads_path, 
                                      fill=TRUE, 
                                      sep = "\t",
                                      header=FALSE) %>% as_tibble() %>% 
  setNames( c("classified", "seq_id", "tax_id", "score", "taxid_score", "taxon_names" ))

if(nrow(df_classed_reads)==0) {
  write_tsv(tibble(), outfile)
  quit(save="no", status=0)
}

df_classed_reads %<>% separate("taxon_names",
                               into=c("superkingdom", "phylum", "order",
                                      "family", "genus", "species", "other"),
                               extra="merge",
                               sep=";") %>% select(-other) %>% 
  mutate_if(is.character, str_trim) %>% 
  mutate_all(na_if,"NA") %>% mutate_all(na_if,"") %>% 
  mutate(species_lower=tolower(str_replace_all(species, " ", "_")) )

df_coverage <- data.table::fread(file = coverage_path, 
                                 fill=TRUE, 
                                 sep = "\t",
                                 header=TRUE) %>% as_tibble()


df_coverage %<>% separate("#ID", into = c("seq_id", "flag", "multi", "len"),
                          extra = "merge", sep = " ") 

df_reads_and_cov <- df_classed_reads %>% 
  left_join(df_coverage, by="seq_id") %>% 
  select(c(seq_id, tax_id, superkingdom, phylum, order, family, genus, species,
           species_lower, Plus_reads, Minus_reads ) )

df_readcount <- 
  df_reads_and_cov %>% 
  rowwise() %>%
  mutate(read_count = sum(Plus_reads, Minus_reads, na.rm = TRUE)) %>% 
  mutate(read_count=ifelse(read_count==0, 1, read_count)) 

df_krona <- 
  df_readcount %>% 
  group_by(superkingdom, phylum, order, family, genus, species) %>% 
  tally() %>% 
  select(n, everything()) %>% 
  write_tsv(outfile,col_names = FALSE)



