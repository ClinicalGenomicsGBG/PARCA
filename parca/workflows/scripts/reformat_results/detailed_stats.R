# Author: Pernilla Ericsson (pernilla.ericsson@gu.se, clinicalgenomics@gu.se)
# Date: 2021-01-07

suppressPackageStartupMessages({
  library(tidyverse);library(magrittr);library(data.table)
} )

read_stats <- function(input_file,
                       classified_reads_name="classified_reads_mincount"){
  df <-data.table::fread(file = input_file, fill=TRUE, sep = "\t",
                         header=TRUE) %>% as_tibble()

  file_basename <- basename(input_file)
  file_basename %<>% str_remove("\\..+")

  if (file_basename == classified_reads_name){
    df %<>% select(superkingdom, seq_id) %>%
      group_by(superkingdom) %>% summarise(count=n(), .groups = 'drop') %>%
      dplyr::rename(type="superkingdom")
  }

  df %<>% mutate(processing_step=file_basename) %>%
    relocate(processing_step, .before = 1)

  return(df)
}

case_control <- snakemake@params[['case_control']]
detailed_stats_out <- snakemake@output[['detailed_stats_out']]

stats_case <- snakemake@input[['stats_case']]
rawpath_case <- snakemake@input[['rawpath_case']]
table_case <- snakemake@input[['table_case']]


stats_list_case <- c(stats_case, table_case)

stats_df_case <- stats_list_case %>% map_dfr(~read_stats(.x)) 

detailed_stats <- stats_df_case

case_name <- table_case %>% str_extract("snakemake_results.+/.{2}_.{3}")
sample_info <- read_tsv(rawpath_case, col_names=TRUE) %>% mutate(processing_step="case_sample_rawpath") %>% bind_rows(tibble(type=case_name, processing_step="case_sample")) 

if (case_control == TRUE){
  
  stats_control <- snakemake@input[['stats_control']]
  rawpath_control <- snakemake@input[['rawpath_control']]
  table_control <- snakemake@input[['table_control']]
  stats_list_control <- c(stats_control, table_control)
  
  stats_df_control <- stats_list_control %>% map_dfr(~read_stats(.x))
  
  detailed_stats <- 
    full_join(stats_df_case, stats_df_control, by=c("processing_step", "type"), 
            suffix = c("_case","_control") ) %>% 
    replace_na(list(type="NA", count_case=0, count_control=0)) 
  
  control_name <- table_control %>% str_extract("snakemake_results.+/.{2}_.{3}")
  sample_info_control <- read_tsv(rawpath_control, col_names=TRUE) %>% mutate(processing_step="control_sample_rawpath") %>% bind_rows(tibble(type=control_name, processing_step="control_sample")) 
  
  sample_info %<>% bind_rows(sample_info_control)
  
}

detailed_stats <- 
  sample_info %>% 
  bind_rows(detailed_stats) %>% 
  mutate(across(!matches("count"), ~as.character(.) )) %>% 
  mutate(across(!matches("count"), ~replace_na(., "NA") )) %>% 
  dplyr::relocate(processing_step, type, .before=1)

write_tsv(detailed_stats, detailed_stats_out, col_names = TRUE)


# case_control <- TRUE
# detailed_stats_out <- "/Users/pernillaericsson/Desktop/parca_tests/detailed_stats.tsv"
# 
# stats_case <- c(
#   raw_reads_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage1/samples/count_raw_reads.txt",
#   trimmed_reads_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage1/trimming/count_bbduk_trimmed_reads.txt",
#   kmer_input_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage2/kmer_input/count_kmer_input.txt",
#   kraken_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage3/kraken/count_kraken_filtered_classified.txt",
#   kaiju_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage3/kaiju/count_kaiju_filtered_classified.txt",
#   species_genus_higher_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage4/count_species_genus_higher.txt",
#   kmer_doublets_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage4/count_doublets.txt",
#   kmer_singletons_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage4/count_singletons.txt",
#   detected_and_missing_slices_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage5/count_detected_missing_slices.txt",
#   subset_blast_reads_taxids_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage6/count_reads_taxid_SubsetBLAST.txt",
#   kmer_subset_blast_reads_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage6/count_kmer_SubsetBLAST.txt",
#   nt_blast_reads_taxids_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage7/count_reads_taxid_BLASTnt.txt",
#   kmer_subset_blast_nt_blast_reads_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage7/count_kmer_SubsetBLAST_BLASTnt.txt",
#   classified_reads_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/PE_RNA/stage8/tableview/classified_reads_mincount.tsv")
# 
# rawpath_case <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage1/rawpath.tsv"
# 
# table_case <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/PE_RNA/stage8/tableview/classified_reads_mincount.tsv"


# 
# stats_control <- c(
#   raw_reads_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage1/samples/count_raw_reads.txt",
#   trimmed_reads_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage1/trimming/count_bbduk_trimmed_reads.txt",
#   kmer_input_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage2/kmer_input/count_kmer_input.txt",
#   kraken_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage3/kraken/count_kraken_filtered_classified.txt",
#   kaiju_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage3/kaiju/count_kaiju_filtered_classified.txt",
#   species_genus_higher_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage4/count_species_genus_higher.txt",
#   kmer_doublets_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage4/count_doublets.txt",
#   kmer_singletons_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage4/count_singletons.txt",
#   detected_and_missing_slices_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage5/count_detected_missing_slices.txt",
#   subset_blast_reads_taxids_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage6/count_reads_taxid_SubsetBLAST.txt",
#   kmer_subset_blast_reads_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage6/count_kmer_SubsetBLAST.txt",
#   nt_blast_reads_taxids_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage7/count_reads_taxid_BLASTnt.txt",
#   kmer_subset_blast_nt_blast_reads_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage7/count_kmer_SubsetBLAST_BLASTnt.txt",
#   classified_reads_case = "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/PE_RNA/stage8/tableview/classified_reads_mincount.tsv")
# 
# rawpath_control <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage1/rawpath.tsv"
# 
# table_control <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/PE_RNA/stage8/tableview/classified_reads_mincount.tsv"
