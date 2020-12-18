# Author: Pernilla Ericsson (pernilla.ericsson@gu.se, clinicalgenomics@gu.se)
# Date: 2020-12-14

suppressPackageStartupMessages({
  library(tidyverse);library(magrittr);library(data.table)
} )

# classed_reads_path <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/200819_demo/snakemake_results_a/SE_RNA/stage8/all_classed_read_taxid_names.txt"
# coverage_path <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/200819_demo/snakemake_results_a/SE_RNA/stage2/pileup/bbmap_cov.txt"
# mincount <- 9
# outfile_readcount <- "/Users/pernillaericsson/Desktop/readcount.tsv"
# outfile_krona <- "/Users/pernillaericsson/Desktop/readcount_krona.tsv"

classed_reads_path <- snakemake@input[['all_classed_read_taxid_names']]
coverage_path <- snakemake@input[['cov']]

mincount <- snakemake@params[['mincount']]

outfile_readcount <- snakemake@output[['readcount']]
outfile_krona <- snakemake@output[['readcount_krona']]

df_classed_reads <- data.table::fread(file = classed_reads_path, 
                                      fill=TRUE, 
                                      sep = "\t",
                                      header=FALSE) %>% as_tibble() %>% 
  setNames( c("classified", "seq_id", "tax_id", "score", "taxid_score", "taxon_names" ))

if(nrow(df_classed_reads)==0) {
  write_tsv(tibble(), outfile_readcount)
  write_tsv(tibble(), outfile_krona)
  quit(save="no", status=0)
}

df_classed_reads %<>% separate("taxon_names",
                               into=c("superkingdom", "phylum", "order",
                                      "family", "genus", "species", "other"),
                               extra="merge",
                               sep=";") %>% select(-other) %>% 
  mutate_if(is.character, str_trim) %>% 
  mutate_all(na_if,"NA") %>% mutate_all(na_if,"") %>% 
  mutate(sgft = pmap_chr(.l = list(species, genus, family, as.character(tax_id) ),
                               .f = ~ first(discard(c(s=..1, g=..2 , f=..3, t=..4), ~is.na(.x) ) ) ) ) %>% 
  mutate(sgft_lower=tolower(str_replace_all(sgft, " ", "_")) ) %>% 
  mutate(tax_id_sgft_lower=paste(tax_id,sgft_lower,sep="_") )


df_coverage <- data.table::fread(file = coverage_path, 
                                 fill=TRUE, 
                                 sep = "\t",
                                 header=TRUE) %>% as_tibble()


df_coverage %<>% separate("#ID", into = c("seq_id", "flag", "multi", "len"),
                          extra = "merge", sep = " ") 

df_reads_and_cov <- 
  df_classed_reads %>% 
  left_join(df_coverage, by="seq_id") 

df_readcount <- 
  df_reads_and_cov %>%
  select(c( classified, seq_id, tax_id, score, taxid_score, superkingdom, 
            phylum, order, family, genus, species, sgft, sgft_lower,
            tax_id_sgft_lower, Plus_reads, Minus_reads ) ) %>% 
  rowwise() %>%
  mutate(read_count = sum(Plus_reads, Minus_reads, na.rm = TRUE)) %>% 
  mutate(read_count=ifelse(read_count==0, 1, read_count)) %>% ungroup() %>% 
  select(-c(Plus_reads, Minus_reads)) %>% 
  write_tsv(outfile_readcount,col_names = TRUE)

df_krona <- 
  df_readcount %>% 
  group_by(superkingdom, phylum, order, family, genus, species) %>% 
  summarize(readcount=sum(read_count)) %>% 
  select(readcount, everything()) %>% 
  arrange(desc(readcount)) %>% 
  filter(readcount > mincount) %>% 
  write_tsv(outfile_krona,col_names = FALSE)



