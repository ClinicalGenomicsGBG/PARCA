


suppressPackageStartupMessages({
  library(tidyverse);library(magrittr); library(ShortRead) #library(data.table)
} )

# trimmed_reads_merged <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1/snakemake_results_sample_1/PE_RNA/stage1/trimming/unmerged_reads.fastq"
# trimmed_reads_unmerged <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1/snakemake_results_sample_1/PE_RNA/stage1/trimming/merged_reads.fastq"
# 
# organism_tableview <- "/Users/pernillaericsson/Desktop/parca_tests/organism_dir/taxid_7962_kingdom_Eukaryota.tsv" # "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1/snakemake_results_sample_1/PE_RNA/stage8/readcount.tsv"
# 
# SE_or_PE <- "PE"
# negate_query <- "FALSE"
# 
# fastq_out <- "/Users/pernillaericsson/Desktop/parca_tests/taxid_test.fastq"

trimmed_reads <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_2/SE_DNA/stage1/trimming/trimmed_reads.fq"
organism_tableview <- "/Users/pernillaericsson/Documents/medair1/medstore/logs/pipeline_logfiles/parca/20201202_run_1_v2_all/snakemake_results_sample_2/SE_DNA/stage8/tableview/organism_dir/organism_100272.tsv"


organism_tableview <- snakemake@input[['organism_tableview']]

SE_or_PE <- snakemake@params[['SE_or_PE']]
negate_query <- snakemake@params[['negate_query']]

fastq_out <- snakemake@output[['fastq_out']]


# Read in the read count stats with regard to PE or SE
if (SE_or_PE == "SE") {
  trimmed_reads <- snakemake@input[['trimmed_reads']]
  fastqs <- c(trimmed_reads)
} else if (SE_or_PE == "PE"){
  trimmed_reads_unmerged <- snakemake@input[['trimmed_reads_unmerged']]
  trimmed_reads_merged <- snakemake@input[['trimmed_reads_merged']]
  fastqs <- c(trimmed_reads_unmerged, trimmed_reads_merged)
} else {
  quit(save="no", status=0)
}


append_to_fastq <- function(fastq_path, ids, outputFile, negate_choice) {
  # Function to append reads to a fastq file.
  fastq <- ShortRead::readFastq(fastq_path)
  i_query <- ids %>% map_chr(~{glue::glue("\\b{.x}\\b")}) %>% paste(collapse = "|")
  subset <- fastq[ str_detect(ShortRead::id(fastq), i_query, negate = negate_choice) ]
  ShortRead::writeFastq(object = subset, file = outputFile, mode="a", compress = FALSE)
}

classified_df <-data.table::fread(file = organism_tableview, fill=TRUE, sep = "\t",
                                  header=TRUE) %>% as_tibble()

read_ids <- classified_df %>% pull(seq_id)

fastqs %>% map(~append_to_fastq(fastq_path=.x, ids=read_ids, outputFile = fastq_out, negate_choice = negate_query))
