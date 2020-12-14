library(tidyverse)
library(magrittr)
library(data.table)

classed_reads_path <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/200819_demo/snakemake_results_a/SE_RNA/stage8/all_classed_read_taxid_names.txt"

coverage_path <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/200819_demo/snakemake_results_a/SE_RNA/stage2/pileup/bbmap_cov.txt"

df_classed_reads <- data.table::fread(file = classed_reads_path, 
                            fill=TRUE, 
                            sep = "\t",
                            header=FALSE) %>% 
  setNames(c("classified", "seq_id", "tax_id", "score", "taxid_score", "taxon_names" )) %>% 
  as_tibble()

df_coverage <- data.table::fread(file = coverage_path, 
                                      fill=TRUE, 
                                      sep = "\t",
                                      header=TRUE) %>%
  as_tibble()

df_coverage %<>% separate("#ID", into = c("seq_id", "flag", "multi", "len"), extra = "merge", sep = " ") 

df_reads_and_cov <- df_classed_reads %>% 
  left_join(df_coverage, by="seq_id") %>% 
  select(classified, seq_id, tax_id, taxon_names, Plus_reads, Minus_reads)

df_reads_and_cov %>% 
  rowwise() %>%
  mutate(read_count = sum(Plus_reads, Minus_reads, na.rm = TRUE)) %>% 
  mutate(read_count=ifelse(read_count==0, 1, read_count)) 
  

