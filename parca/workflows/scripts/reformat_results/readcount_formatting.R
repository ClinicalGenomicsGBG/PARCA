# Author: Pernilla Ericsson (pernilla.ericsson@gu.se, clinicalgenomics@gu.se)
# Date: 2020-12-21

suppressPackageStartupMessages({
  library(tidyverse);library(magrittr);library(data.table)
} )

# classed_reads_path <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/200819_demo/snakemake_results_a/SE_RNA/stage8/all_classed_read_taxid_names.txt"
# mincount <- 9
# DNA_or_RNA <- "RNA"
# coverage_path <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/200819_demo/snakemake_results_a/SE_RNA/stage2/pileup/bbmap_cov.txt"
# outfile_readcount <- "/Users/pernillaericsson/Desktop/readcount.tsv"
# outfile_krona <- "/Users/pernillaericsson/Desktop/readcount_krona.tsv"

classed_reads_path <- snakemake@input[['all_classed_read_taxid_names']]

mincount <- snakemake@params[['mincount']]
DNA_or_RNA <- snakemake@params[['DNA_or_RNA']]

outfile_readcount <- snakemake@output[['readcount']]
outfile_krona <- snakemake@output[['readcount_krona']]


df_classed_reads <- data.table::fread(file = classed_reads_path, fill=TRUE, 
                                      sep = "\t", header=FALSE) %>% 
  as_tibble() %>% 
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
  mutate(organism = pmap_chr(.l = list(species, genus, family, as.character(tax_id) ),
                             .f = ~first(discard(c(s=..1, g=..2 , f=..3, t=..4), ~is.na(.x) ) ) ) )

if (DNA_or_RNA == "RNA") {
  coverage_path <- snakemake@input[['cov']]

  df_coverage <- data.table::fread(file = coverage_path, 
                                   fill=TRUE, 
                                   sep = "\t",
                                   header=TRUE) %>% as_tibble()
  df_coverage %<>% separate("#ID", into = c("seq_id", "flag", "multi", "len"), extra = "merge", sep = " ") 
  
  df_reads_and_cov <- 
    df_classed_reads %>% 
    left_join(df_coverage, by="seq_id") 
  
  df_readcount <- 
    df_reads_and_cov %>%
    select(c( classified, seq_id, tax_id, score, taxid_score, superkingdom, 
              phylum, order, family, genus, species, organism, 
              Plus_reads, Minus_reads ) ) %>% 
    mutate(Plus_reads=as.integer(Plus_reads), Minus_reads=as.integer(Minus_reads)) %>% 
    rowwise() %>%
    mutate(read_count = sum(Plus_reads, Minus_reads, na.rm = TRUE)) %>% 
    mutate(read_count=ifelse(read_count==0, 1, read_count)) %>% ungroup() %>% 
    select(-c(Plus_reads, Minus_reads)) %>% 
    write_tsv(outfile_readcount,col_names = TRUE) 
  
  } else if (DNA_or_RNA == "DNA"){
  df_readcount <- 
    df_classed_reads %>%
    mutate(read_count=1) %>% 
    write_tsv(outfile_readcount,col_names = TRUE)
} else {
  quit(save="no", status=0)
}

df_krona <- 
  df_readcount %>% 
  group_by(superkingdom, phylum, order, family, genus, species) %>% 
  summarize(readcount=sum(read_count)) %>% 
  select(readcount, everything()) %>% 
  arrange(desc(readcount)) %>% 
  filter(readcount > mincount) %>% 
  write_tsv(outfile_krona,col_names = FALSE)
