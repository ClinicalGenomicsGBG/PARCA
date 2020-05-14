# Script for moving all reads that were classed below species, up to a taxid corresponding to a species or higher.
# Also an xml file containing the taxids of reads that should be considered as primates is used. If it is empty, nothing will be filtered.
# If either the taxid lineage file or the filtered blast results file is empty then the blast file will be returned with extra column names or an empty file will be returned respectively.
# The new result file (variable species_and_above) will contain the columns "seq_id,	tax_id,	score,	rank,	organism,	type", if the blast file is not empty
# 
# Author: Pernilla Ericsson (pernilla.ericsson@gu.se)
# Date: 2020-05-12

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(data.table) 
  library(methods)
  library(XML)
} )

best_blast <- snakemake@input[['best_blast']]
tax_id_lineage <- snakemake@input[['tax_id_lineage']]

primates_file <- snakemake@params[['primates_file']]
blast_type <- snakemake@params[['blast_type']]

species_and_above <- snakemake@output[['species_and_above']]
count_reads <- snakemake@output[['count_reads_tax_ids']]

# Read the filtered blast file.
df_blast <- data.table::fread(file = best_blast, 
                  fill=TRUE, 
                  sep = "\t",
                  header=TRUE) %>% 
  as_tibble()

# Read the taxid to lineage file
lineage_df_raw <- data.table::fread(file = tax_id_lineage, 
                                    fill=TRUE, 
                                    sep = "\t",
                                    header=FALSE) %>% 
  as_tibble() %>% 
  setNames(c("tax_id","rank", "taxonomy") )

# Check the size of the primates taxid file
primates_file_info = file.info(primates_file, extra_cols = FALSE)$size %>% as.double()

# Quit program if either blast result file is empty or lineage file is empty.
if (nrow(df_blast)==0) {
  write_tsv(tibble(), species_and_above)
  write_tsv(tibble(type="Total", count=0), count_reads)
  quit(save = "no", status = 0)
} else if (nrow(lineage_df_raw)==0) {
  df_blast %<>% 
    dplyr::rename("seq_id"=qseqid, "tax_id"=staxids) %>% 
    mutate(rank="s", organism="SPECIES", type=blast_type) %>% 
    select(seq_id,	tax_id,	score,	rank,	organism,	type)
  write_tsv(df_blast, species_and_above)
  write_tsv(tibble(type="Total", count=nrow(df_blast)), count_reads)
  quit(save = "no", status = 0)
}

# If primates file is empty, then return an empty list. If it is not empty, read the xml file and return a list of taxids.
if (primates_file_info==0) {
  primates_list <- c()
} else{
  xml_result <- xmlParse(file = primates_file)
  df_ids <- xmlToDataFrame(nodes = xmlChildren(xmlRoot(xml_result)[["IdList"]])) 
  primates_list <- df_ids %>% pull(text) %>% as.character()
}

# Reformat the lineage file.
lineage_df <- 
  lineage_df_raw %>% 
  mutate(tax_id=as.character(tax_id)) %>% 
  separate(col = taxonomy,into = c("kingdom", "phylum", "class", "order", "family", "genus","species"), sep = ";") %>% 
  mutate_all(na_if,"")

# Reformat the blast results file.
df_blast %<>% 
  dplyr::rename("seq_id"=qseqid, "tax_id"=staxids) %>% 
  mutate(tax_id=as.character(tax_id))

# Add lineage to blast results.
df_blast_lineage <- left_join(df_blast, lineage_df, by="tax_id")

# Remove unneccessary variables to save space.
rm(df_blast,lineage_df_raw,lineage_df,df_ids)

# Move up all reads classed to a taxid lower than species, to species or higher.
df_blast_lineage %<>% 
  replace_na(list(rank = "no rank")) %>% 
  mutate(final_tax_id= ifelse(
    !is.na(species),
    species,
    ifelse(!str_detect(rank, "species"),
           tax_id,
           paste(tax_id, kingdom, phylum, class, order, family, genus,species, sep="__") %>%
             str_remove_all(pattern = "NA__|__NA|NA") %>% str_remove("\\w+__") 
            ) )
  )  %>% 
  select(seq_id, final_tax_id, score)

NAs_in_final <- 
  df_blast_lineage %>% pull(final_tax_id) %>% is.na() %>% sum()

# If NAs are found in the output file there must be a bug somewhere...
if (NAs_in_final != 0){
  print("Found NA in final taxid!!")
  quit(save = "no", status = 1)
}

# Reformat and filter the output file to change all taxids in primates file to 9606.
df_blast_lineage_clean <- 
  df_blast_lineage %>% 
  dplyr::rename("tax_id"=final_tax_id) %>% 
  mutate(rank = "s", organism="SPECIES", type=blast_type) %>% 
  select(seq_id,	tax_id,	score,	rank,	organism,	type) %>% 
  filter(!tax_id %in% "0") %>% 
  mutate(tax_id=ifelse(tax_id %in% "10870", "1977402", tax_id),
         tax_id=ifelse(tax_id %in% primates_list, "9606", tax_id)
         )

write_tsv(df_blast_lineage_clean, species_and_above, col_names = TRUE)

# Print some stats! 
no_tax_ids <- df_blast_lineage_clean %>% pull(tax_id) %>% unique() %>% length()
no_reads <- nrow(df_blast_lineage_clean)

stats <- tibble(type=c("tax ids","reads"), count=c(no_tax_ids, no_reads))
write_tsv(stats, count_reads, col_names = TRUE)


# best_blast <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage6/best_blast.txt"
# tax_id_lineage <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage6/best_blast_tax_id_lineage.txt"
# blast_type <- "SubsetBLAST"
# primates_file <- "/Users/pernillaericsson/Documents/medair1/medstore/Development/Metagenomics/PathFinder43b/databases/primates_taxids.txt"
# species_and_above <- "/Users/pernillaericsson/Desktop/slice_blasted_reads_species_and_above"
# count_reads <- "/Users/pernillaericsson/Desktop/count_slice_blasted_reads"
# empty_file <- "/Users/pernillaericsson/Desktop/empty_test.txt"
# 
# 
# # df_blast_lineage2 <-
#   df_blast_lineage %>%
#   bind_rows(structure(list(seq_id = "LH54F:00015:06573", tax_id = "296587",
#                            score = 197.4, rank = "species", kingdom = "2759", phylum = "3041",
#                            class = "1035538", order = NA, family =NA, genus = NA,
#                            species = NA), row.names = c(NA, -1L), class = c("tbl_df",
#                                                                                   "tbl", "data.frame")),.)