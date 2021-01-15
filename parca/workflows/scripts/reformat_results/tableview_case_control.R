suppressPackageStartupMessages({
  library(tidyverse);library(magrittr);library(data.table);library(ShortRead)
} )

# tableview_case <- "/Users/pernillaericsson/Desktop/parca_tests/tableview.tsv"
# tableview_control<- "/Users/pernillaericsson/Desktop/parca_tests/tableview.tsv"
# 
# tableview_case_control_out <- "/Users/pernillaericsson/Desktop/parca_tests/tableview_case_control"

tableview_case <- snakemake@input[['case']]
tableview_control <- snakemake@input[['control']]
tableview_case_control_out  <- snakemake@output[['case_control']]


df_case <- data.table::fread(input = tableview_case, fill = TRUE, sep="\t", header=TRUE) %>% 
  as_tibble() %>% mutate(tax_id=as.character(tax_id))

df_control <- data.table::fread(input = tableview_control, fill = TRUE, sep="\t", header=TRUE) %>% 
  as_tibble() %>% mutate(tax_id=as.character(tax_id))

df_case %>% full_join(df_control,
                      by=c("superkingdom", "organism", "tax_id"),
                      suffix = c("_case","_control")) %>% 
  mutate(across(where(is.numeric), ~replace_na(., 0) )) %>% 
  mutate(across(where(is.character), ~replace_na(., "NA") )) %>% 
  write_tsv(tableview_case_control_out, col_names = TRUE)
