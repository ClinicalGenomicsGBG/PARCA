#
#

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(data.table) 
} )

# all_downloaded <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage5/downloadblastslices/downloaded_slices.txt"
# splitaccdump_dir <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/ref/splitaccdump"
# 
# created_slice_dir <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/ref/"

all_downloaded <-snakemake@input[['all_downloaded']]

splitaccdump_dir <- snakemake@params[['splitaccdump_dir']]

created_slice_dir <- snakemake@output[['created_slice_dir']]


# Read file with the downloaded taxas. 
df_downloaded <-  data.table::fread(file = all_downloaded) %>% 
  mutate_all(as.character)

if(nrow(df_downloaded)==0){
  system(paste("if [ ! -d", created_slice_dir," ]; then mkdir", created_slice_dir, ";fi;") )
  quit(save = "no", status = 0)
}

colnames(df_downloaded) <- c("file", "taxid") 

# Get all dmp files to be filtered.
acc_dmp_files <- list.files(path=splitaccdump_dir, pattern = "*.dmp$")

# Function to extract downloaded taxas from dmp file with accessions.
find_accessions_to_taxids <- function(dmp_file) {
  dmp_file_name <- paste0(splitaccdump_dir,"/",dmp_file)
  df_dmp <- data.table::fread(file = dmp_file_name,header=FALSE) %>% as_tibble()
  colnames(df_dmp) <- c("accession","accession.version", "taxid", "gi")
  df_dmp %<>% dplyr::select(accession.version,taxid) %>% 
    mutate_all(as.character)
  df_tmp <- inner_join(df_downloaded, df_dmp,by="taxid")
}

# Bind all df_tmp (dmp files filtered for downloaded taxas)
df_taxid_acc <- 
  acc_dmp_files %>% 
  map_dfr(~find_accessions_to_taxids(.x))


downloaded_taxids <- df_downloaded %>% pull(file) %>% unique()

# Function to save the df with the name of the file (main taxa)
save_by_taxid <- function(current_taxid){
  out_name <- paste0(created_slice_dir,"/",current_taxid)
  filtered_slice <- df_taxid_acc %>% 
    filter(file == current_taxid) %>% 
    dplyr::select(accession.version)
  write_tsv(filtered_slice, out_name,col_names = FALSE) 
  }


system(paste("if [ ! -d",created_slice_dir," ]; then \nmkdir",created_slice_dir, ";fi;") )
downloaded_taxids %>% 
  map(~save_by_taxid(.x))
