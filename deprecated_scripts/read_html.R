#!/usr/bin/env Rscript
library(tidyr)
library(purrr)
library(openxlsx)
require(rvest)
require(magrittr)
require(tidyverse)

html_reader <- function(html_file, parsing=TRUE) {
  # require(rvest)
  # require(magrittr)
  # require(tidyverse)
  
  webpage <- read_html(html_file)
  
  tbls_ls <- webpage %>%
    html_nodes("table") %>% 
    html_table(fill = TRUE) 
  
  if(parsing==TRUE){
  res <-
    tbls_ls %>% map(~{ 
      df <- .x 
      type_name <- df %>% colnames() %>%  magrittr::extract(1) 
      colnames(df) <- paste("col", 1:ncol(df),sep = "_") 
      c_names <- df %>% slice(1) %>% unname() 
      df %<>% slice(-1)
      colnames(df) <- c_names
      df %<>% add_column(Type= type_name,.before = 1)
    }) %>% bind_rows() %>% as_tibble()
  
  return(res)
  } else {
    return(tbls_ls)
  }
  
}

create_sample_summary <- function(basePath="/seqstore/webfolders/PathFinder", sampleFolder, tableFile="tableview.html", statsFile="stats.html"){
  samples <- sampleFolder %>% str_extract("[GMS]-\\d+_[GMS]-\\d+|[GMS]-\\d+")
  print(samples)
  tableFilePath <- paste0(basePath,"/", sampleFolder,"/", tableFile)
  table <- html_reader(tableFilePath)
  df <- 
    table %>% 
    mutate(Hits=str_replace_all(Hits,"[:punct:]", "")) %>% 
    mutate(Hits=str_replace_all(Hits," ", "_"))  %>% 
    mutate(type__taxid__organism=paste0(Type,"__",Taxid,"__",Hits)) %>% 
    mutate(type__taxid__organism, type__taxid__organism=tolower(type__taxid__organism)) 
  
  
  statsFilePath <- paste0(basePath,"/", sampleFolder,"/", statsFile)
  stats <- html_reader(statsFilePath, parsing = FALSE)
  trimmed_stats <- stats[[2]]
  colnames(trimmed_stats) <- c("sample_type","trimmed_reads")
  trimmed_stats %<>% slice(-1)
  
  rna_reads <- trimmed_stats %>% filter(sample_type=="RNA:") %>% pull(trimmed_reads)
  nc <- trimmed_stats %>% filter(sample_type=="NC_RNA:") %>% nrow()
  
  if (nc==0) {
    rna_df <- 
      df %>% 
      select(type__taxid__organism, "Reads (RNA)")
    
    rna_sample <-  samples
    colnames(rna_df) <- c("type__taxid__organism",rna_sample)
  
    reads_row <- c("reads_after_trimming",rna_reads)
    
  } else{
    rna_df <- 
      df %>% 
      select(type__taxid__organism, "Reads (RNA)", "Reads (NC_RNA)") 
    
    rna_sample <- samples %>% str_extract("[GMS]-\\d+(?=_)")
    nc_rna_sample <- samples %>% str_extract("(?<=_)[GMS]-\\d+")
    colnames(rna_df) <- c("type__taxid__organism",rna_sample,nc_rna_sample)
    
    nc_rna_reads <- trimmed_stats %>% filter(sample_type=="NC_RNA:") %>% pull(trimmed_reads)
    reads_row <- c("reads_after_trimming",rna_reads, nc_rna_reads)
    
  }
  
  rna_df <- rbind(reads_row, rna_df)
  return(rna_df)
  
}


# base="/seqstore/webfolders/PathFinder"
# df_full_out="/Users/pernillaericsson/Desktop/df_full.xlsx"
# df_full_sampleRows_out="/Users/pernillaericsson/Desktop/df_full_sampleRows.xlsx"

#directories <- list.dirs(path=base,full.names = FALSE)
# directories %<>% 
#   na_if("") %>%
#   na.omit() %>% 
#   c()

#directories <- c("2020-01-29_G-261_G-23_15_22_8367", "2020-01-31_G-393_01_31_1769" )


args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file of directories).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[1] = ""
  args[2] = "/seqstore/webfolders/PathFinder"
  args[3] = "df_full.xlsx"
  args[4] = "df_full_sampleRows.xlsx"
  
}
## program...
directories_file <- args[1] 
folder <- args[2] 
df_full_out <- args[3] 
df_full_sampleRows_out <- args[4] 

#directories_file= "/Users/pernillaericsson/Documents/medair1/seqstore/webfolders/PathFinder/dir.txt"
directories <- read_csv(directories_file,col_names = "dir") %>% pull(dir)

df_full_all_list <- directories %>%
  map(~create_sample_summary(basePath=folder, sampleFolder = .x))

df_full_all <- 
  df_full_all_list %>% 
  reduce(full_join, by = "type__taxid__organism")

df_full_nona <-  df_full_all %>% mutate_all(~(replace_na(.x,0)))

full_nona_sampleRows <- df_full_nona %>% t()
c_names = full_nona_sampleRows[1,]
colnames(full_nona_sampleRows) <- c_names
df_full_nona_sampleRows <- full_nona_sampleRows[-1,]
df_full_nona_sampleRows

print(dim(df_full_nona_sampleRows)) #144 11227

write.xlsx(df_full_nona,file=df_full_out,colNames=TRUE, rowNames=FALSE)
write.xlsx(df_full_nona_sampleRows,file=df_full_sampleRows_out,colNames=TRUE, rowNames=TRUE)



case_control_file <- "/Users/pernillaericsson/Documents/medair1/home/xerpey/pathfinder_runs/140_samples/controls_cases.xlsx"
case_control <- read.xlsx(case_control_file)
case_controls_long <- case_control %>% pivot_longer(everything(), names_to = "type", values_to = "sample") %>% na.omit()

case_controls_long %<>%
  mutate(sample=paste0(str_extract(sample,"[GMS]"), "-", str_extract(sample,"\\d+")) ) %>% 
  mutate(type=str_replace(type, "Controls", "control")) %>% 
  mutate(type=str_replace(type, "Cases", "case"))

full <- read.xlsx("/Users/pernillaericsson/Documents/medair1/home/xerpey/pathfinder_runs/140_samples/df_full/full_sampleRows.xlsx")
colnames(full)[1] <- "sample"
dim(full)

full_casecontrol <- left_join(full,case_controls_long,by="sample")
dim(full_casecontrol)

full_casecontrol_reordered <- full_casecontrol %>% select(sample, type, reads_after_trimming,everything())
dim(full_casecontrol_reordered)

write.xlsx(full_casecontrol_reordered,file="~/Desktop/detected_organisms_parca.xlsx",colNames=TRUE, rowNames=FALSE)