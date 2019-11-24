
library(tidyverse);library(magrittr);library(openxlsx)

#specification_file ="~/Documents/projects/pathfinder/controls_cases.xlsx"
#all_samples_path="~/Documents/projects/pathfinder/samples"
#index_name="index_140_research_samples"
#illumina="yes"
#outfile_name="~/Desktop/output.txt"

####
specification_file ="~/pathfinder_runs/140_samples/controls_cases.xlsx"
all_samples_path="~/pathfinder_runs/140_samples/samples"
index_name="~/pathfinder_runs/140_samples/index_140_research_samples.html"
illumina="yes"
outfile_name="~/pathfinder_runs/140_samples/snakemake_interleave_files.txt"
###


specification <- read.xlsx(specification_file)
find_file_paths <- function(sample_id, file_path) {
  id <- str_extract(sample_id, "^[:alpha:]")
  number <- str_extract(sample_id, "\\d+")
  #print(number)
  
  sample_pattern <- paste0("^",id,"[_-]{0,1}",number,"_")
  #print(sample_pattern)
  
  if(is.na(number[1])){
    sample_path <- NA
    print(sample_path)
  } else {
    sample_path <- list.files(path=file_path, pattern=sample_pattern) 
    #print(sample_path)
  }
  
  
  return(sample_path[1])
}


launch_specification <- specification %>%
  mutate(Control_path = map_chr(Controls,~find_file_paths(.x,all_samples_path))) %>% 
  mutate(Cases_path = map_chr(Cases,~find_file_paths(.x,all_samples_path))) 

head(launch_specification)


id_specification <- launch_specification %>% 
  separate(Control_path, c("control_id", "control_rest"), "_(R1|R2)_val_1.fq", extra = "merge") %>% 
  separate(Cases_path, c("case_id", "case_rest"), "_(R1|R2)_val_1.fq", extra = "merge")

control_ids <- id_specification %>% pull(control_id) %>% na.omit(.)
case_ids <-id_specification %>% pull(case_id) %>% na.omit(.)


config_ids <- c(control_ids, case_ids)

snakemake_ids <- config_ids %>% paste("-",., collapse = "\n")

outfile<-file(outfile_name)
writeLines(snakemake_ids, outfile)
close(outfile)



