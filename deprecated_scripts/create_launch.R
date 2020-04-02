
library(tidyverse);library(magrittr)
library(openxlsx)

#specification_file ="~/Documents/projects/pathfinder/controls_cases.xlsx"
#all_samples_path="~/Documents/projects/pathfinder/samples"
#index_name="index_140_research_samples"
#illumina="yes"
#outfile_name="~/Desktop/output.txt"


#specification_file ="~/pathfinder_runs/140_samples/controls_cases.xlsx"
specification_file = "/apps/bio/repos/pathfinder_python/references/controls_cases.xlsx"
all_samples_path="/medstore/Development/Metagenomics/projects/140_cervix_samples_pathfinder/interleaved_samples"
index_name="index_cervix_samples_v3.html"
illumina="yes"
#outfile_name="~/pathfinder_runs/140_samples/runvalidate_140_research_samples2.pl"
outfile_name = "/apps/bio/repos/pathfinder_python/runvalidate_140_research_samples_v3.pl"

specification <- read.xlsx(specification_file)

find_file_paths <- function(sample_id, file_path) {
  id <- str_extract(sample_id, "^[:alpha:]")
  number <- str_extract(sample_id, "\\d+")
  
  sample_pattern <- paste0("^",id,"[_-]{0,1}",number,"_")
  
  if(is.na(number[1])){
    sample_path <- NA
    print(sample_path)
  } else {
    sample_path <- list.files(path=file_path, pattern=sample_pattern, full.names = TRUE) 
  }
  
  
  return(sample_path[1])
}


launch_specification <- specification %>%
  mutate(Control_path = map_chr(Controls,~find_file_paths(.x,all_samples_path))) %>% 
  mutate(Cases_path = map_chr(Cases,~find_file_paths(.x,all_samples_path)))

settings <- launch_specification %>% 
  mutate(
         Cases_path = ifelse(is.na(Cases_path),NA, paste0("RNA=",Cases_path) ),
         Control_path = map2_chr(Control_path,Cases_path,~{ifelse(is.na(.y), paste0("RNA=",.x),paste0("NC_RNA=",.x))}),
         Control_path = na_if(Control_path,"NC_RNA=NA"),
         sample_out  = imap_chr(Controls,~{
           
           case <- Cases[.y]
           control <- .x
           control_path <- Control_path[.y]
           case_path <- Cases_path[.y]
           
           case_when(!is.na(control_path) & !is.na(case_path) ~ paste0("SAMPLE=",case,"_",control),
                     is.na(control_path) & !is.na(case_path) ~ paste0("SAMPLE=", case),
                     !is.na(control_path) & is.na(case_path) ~ paste0("SAMPLE=",control)
           )
         })
  ) 

settings %>% 
  select(-c(Controls,Cases)) %>% 
  map_df(~{replace_na(.x,"")}) %>% 
  transmute(launch = paste(sample_out,Cases_path, Control_path, sep=" ")) %>%
  map_df(~{.x %>% as.character() %>% str_trim(side = "both")}) %>%
  mutate(launch = paste0("\'./launch_v2.sh ", launch, " ILLUMINA=", illumina, " WEBPATH=",index_name,"\'")) %>% 
  pull(launch) %>% paste(collapse = ",\n") %>% 
  paste("@array=(",.,");\n\n foreach (@array){\n$time=localtime();\n print \"$time\";\n print \"$_\";\n system(\"$_\");\n}",collapse = "",sep="") %>% 
  write_lines(path = outfile_name,sep = "\n")


