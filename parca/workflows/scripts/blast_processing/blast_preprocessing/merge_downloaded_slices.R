# Script for creating a file with all downloaded taxids with their accessions.
# Author: Pernilla Ericsson (pernilla.ericsson@gu.se)
# Date: 2020-05

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(methods)
  library(XML)
} )

#downloaded_dir <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage5/downloadblastslices/downloaded_slices/" 

downloaded_dir <- snakemake@input[['downloaded_dir']]
all_downloaded <- snakemake@output[['all_downloaded']]

taxonomy_files <- list.files(downloaded_dir, pattern = "\\d+")

read_taxonomy_files <- function(tax_id){
  tax_file <- paste0(downloaded_dir,"/",tax_id)
  
  xml_result <- xmlParse(file = tax_file)
  df_ids <- xmlToDataFrame(nodes = xmlChildren(xmlRoot(xml_result)[["IdList"]]))
  
  df_ids %<>% 
    dplyr::rename("id"=text) %>% 
    mutate(file=tax_id) %>% 
    select(file, id) %>% 
    mutate_all(as.character)
  
  # file_name_in_df <- df_ids %>% filter(id %in% tax_id) %>% nrow()
  # if (file_name_in_df==0) {
  #   df_ids %<>% bind_rows(tibble(file=tax_id, id=tax_id))
  # }
  
  return(df_ids)
}

df_all <- 
  taxonomy_files %>%
  map_dfr(~{read_taxonomy_files(.x)})

df_all %>% write_tsv(path = all_downloaded)