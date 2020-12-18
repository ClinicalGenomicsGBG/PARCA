# Author: Pernilla Ericsson (pernilla.ericsson@gu.se, clinicalgenomics@gu.se)
# Date: 2020-12-14

suppressPackageStartupMessages({
  library(tidyverse);library(magrittr);library(data.table)
} )

# classed_reads_path <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/200819_demo/snakemake_results_a/SE_RNA/stage8/all_classed_read_taxid_names.txt"
# mincount <- 9  
# outfile_readcount <- "/Users/pernillaericsson/Desktop/readcount.tsv"
# outfile_krona <- "/Users/pernillaericsson/Desktop/readcount_krona.tsv"

classed_reads_path <- snakemake@input[['all_classed_read_taxid_names']]

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
  mutate(species_lower=tolower(str_replace_all(species, " ", "_")) )

df_readcount <- 
  df_classed_reads %>%
  mutate(read_count=1) %>% 
  write_tsv(outfile_readcount,col_names = TRUE)

df_krona <- 
  df_readcount %>% 
  group_by(superkingdom, phylum, order, family, genus, species, species_lower) %>% 
  summarize(readcount=sum(read_count)) %>% 
  select(readcount, everything()) %>% 
  arrange(desc(readcount)) %>% 
  filter(readcount > mincount) %>% 
  write_tsv(outfile_krona,col_names = FALSE)



