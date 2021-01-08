# Author: Pernilla Ericsson (pernilla.ericsson@gu.se, clinicalgenomics@gu.se)
# Date: 2021-01-07

suppressPackageStartupMessages({
  library(tidyverse);library(magrittr);library(data.table)
} )

read_stats <- function(input_file,
                       classified_reads_name="classified_reads_mincount"){
  df <-data.table::fread(file = input_file, fill=TRUE, sep = "\t",
                         header=TRUE) %>% as_tibble()
  
  if (nrow(df) == 0){
    return(df)
  }
  
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
# 
# case_control <- TRUE
# detailed_stats_out <- "/Users/pernillaericsson/Desktop/parca_tests/detailed_stats.tsv"
# 
# raw_reads_case <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage1/samples/count_raw_reads.txt"
# trimmed_reads_case <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage1/trimming/count_bbduk_trimmed_reads.txt"
# kmer_input_case <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage2/kmer_input/count_kmer_input.txt"
# kraken_case <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage3/kraken/count_kraken_filtered_classified.txt"
# kaiju_case <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage3/kaiju/count_kaiju_filtered_classified.txt"
# species_genus_higher_case <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage4/count_species_genus_higher.txt"
# kmer_doublets_case <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage4/count_doublets.txt"
# kmer_singletons_case <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage4/count_singletons.txt"
# detected_and_missing_slices_case <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage5/count_detected_missing_slices.txt"
# subset_blast_reads_taxids_case <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage6/count_reads_taxid_SubsetBLAST.txt"
# kmer_subset_blast_reads_case <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage6/count_kmer_SubsetBLAST.txt"
# nt_blast_reads_taxids_case <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage7/count_reads_taxid_BLASTnt.txt"
# kmer_subset_blast_nt_blast_reads_case <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage7/count_kmer_SubsetBLAST_BLASTnt.txt"
# classified_reads_case <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/PE_RNA/stage8/tableview/classified_reads_mincount.tsv"


case_control <- snakemake@params[['case_control']]
detailed_stats_out <- snakemake@output[['detailed_stats_out']]

raw_reads_case <- snakemake@input[['raw_reads_case']]
trimmed_reads_case <- snakemake@input[['trimmed_reads_case']]
kmer_input_case <- snakemake@input[['kmer_input_case']]
kraken_case <- snakemake@input[['kraken_case']]
kaiju_case <- snakemake@input[['kaiju_case']]
species_genus_higher_case <- snakemake@input[['species_genus_higher_case']]
kmer_doublets_case <- snakemake@input[['kmer_doublets_case']]
kmer_singletons_case <- snakemake@input[['kmer_singletons_case']]
detected_and_missing_slices_case <- snakemake@input[['detected_and_missing_slices_case']]
subset_blast_reads_taxids_case <- snakemake@input[['subset_blast_reads_taxids_case']]
kmer_subset_blast_reads_case <- snakemake@input[['kmer_subset_blast_reads_case']]
nt_blast_reads_taxids_case <- snakemake@input[['nt_blast_reads_taxids_case']]
kmer_subset_blast_nt_blast_reads_case <- snakemake@input[['kmer_subset_blast_nt_blast_reads_case']]
classified_reads_case <- snakemake@input[['classified_reads_case']]
  
stats_list_case <- c(raw_reads_case, trimmed_reads_case, kmer_input_case, 
                     kraken_case, kaiju_case, species_genus_higher_case, kmer_doublets_case,
                     kmer_singletons_case, detected_and_missing_slices_case, 
                     subset_blast_reads_taxids_case, kmer_subset_blast_reads_case,
                     nt_blast_reads_taxids_case, kmer_subset_blast_nt_blast_reads_case,
                     classified_reads_case)

stats_df_case <- stats_list_case %>% map_dfr(~read_stats(.x))

detailed_stats <- stats_df_case %>% relocate(processing_step, type, count, .before=1)

case_name <- classified_reads_case %>% str_extract("snakemake_results.+/.{2}_.{3}")
sample_info <- tibble(processing_step=c("case_sample"),
                      type=c(case_name))

if (case_control == TRUE){
  # raw_reads_control <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage1/samples/count_raw_reads.txt"
  # trimmed_reads_control <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage1/trimming/count_bbduk_trimmed_reads.txt"
  # kmer_input_control <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage2/kmer_input/count_kmer_input.txt"
  # kraken_control <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage3/kraken/count_kraken_filtered_classified.txt"
  # kaiju_control <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage3/kaiju/count_kaiju_filtered_classified.txt"
  # species_genus_higher_control <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage4/count_species_genus_higher.txt"
  # kmer_doublets_control <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage4/count_doublets.txt"
  # kmer_singletons_control <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage4/count_singletons.txt"
  # detected_and_missing_slices_control <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage5/count_detected_missing_slices.txt"
  # subset_blast_reads_taxids_control <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage6/count_reads_taxid_SubsetBLAST.txt"
  # kmer_subset_blast_reads_control <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage6/count_kmer_SubsetBLAST.txt"
  # nt_blast_reads_taxids_control <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage7/count_reads_taxid_BLASTnt.txt"
  # kmer_subset_blast_nt_blast_reads_control <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/stats_PE_RNA/stage7/count_kmer_SubsetBLAST_BLASTnt.txt"
  # classified_reads_control <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_1/PE_RNA/stage8/tableview/classified_reads_mincount.tsv"
  
  raw_reads_control <- snakemake@input[['raw_reads_control']]
  trimmed_reads_control <- snakemake@input[['trimmed_reads_control']]
  kmer_input_control <- snakemake@input[['kmer_input_control']]
  kraken_control <- snakemake@input[['kraken_control']]
  kaiju_control <- snakemake@input[['kaiju_control']]
  species_genus_higher_control <- snakemake@input[['species_genus_higher_control']]
  kmer_doublets_control <- snakemake@input[['kmer_doublets_control']]
  kmer_singletons_control <- snakemake@input[['kmer_singletons_control']]
  detected_and_missing_slices_control <- snakemake@input[['detected_and_missing_slices_control']]
  subset_blast_reads_taxids_control <- snakemake@input[['subset_blast_reads_taxids_control']]
  kmer_subset_blast_reads_control <- snakemake@input[['kmer_subset_blast_reads_control']]
  nt_blast_reads_taxids_control <- snakemake@input[['nt_blast_reads_taxids_control']]
  kmer_subset_blast_nt_blast_reads_control <- snakemake@input[['kmer_subset_blast_nt_blast_reads_control']]
  classified_reads_control <- snakemake@input[['classified_reads_control']]

  stats_list_control <- c(raw_reads_control, trimmed_reads_control, kmer_input_control, 
    kraken_control, kaiju_control, species_genus_higher_control, kmer_doublets_control,
    kmer_singletons_control, detected_and_missing_slices_control,
    subset_blast_reads_taxids_control, kmer_subset_blast_reads_control,
    nt_blast_reads_taxids_control, kmer_subset_blast_nt_blast_reads_control,
    classified_reads_control)
  
  stats_df_control <- stats_list_control %>% map_dfr(~read_stats(.x))
  
  detailed_stats <- 
    full_join(stats_df_case, stats_df_control, by=c("processing_step", "type"), 
            suffix = c("_case","_control") ) %>% 
    replace_na(list(type="NA", count_case=0, count_control=0)) %>% 
    relocate(processing_step, type, .before=1)
  
  control_name <- classified_reads_control %>% str_extract("snakemake_results.+/.{2}_.{3}")
  sample_info %<>%  bind_rows(tibble(processing_step=c("control_sample"),
                        type=c(control_name)))

}

detailed_stats <- 
  sample_info %>% 
  bind_rows(detailed_stats)

write_tsv(detailed_stats, detailed_stats_out, col_names = TRUE)
