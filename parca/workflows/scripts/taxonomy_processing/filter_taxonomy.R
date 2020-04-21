#
#Script for setting ranks below genus to genus.
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(data.table) 
  } )

singletons_file <- snakemake@input[['singletons_added_SGF_empty']]
lineage_file <- snakemake@input[['tax_id_lineage']]

singletons_genus_file <- snakemake@output[['singletons_genus']]
singletons_genus_stats_file <- snakemake@output[['read_count_singletons_genus']]

df <- data.table::fread(file = singletons_file, 
                        fill=TRUE, 
                        sep = "\t",
                        header=TRUE) %>% 
  as_tibble() 

df

lineage_df_raw <- data.table::fread(file = lineage_file, 
                             fill=TRUE, 
                             sep = "\t",
                             header=FALSE) %>% 
  as_tibble() %>% 
  setNames(c("tax_id","rank", "taxonomy") )

lineage_df <- 
  lineage_df_raw %>% 
  separate(col = taxonomy,into = c("kingdom", "phylum", "class", "order", "family", "genus","species"), sep = ";") %>% 
  mutate_all(na_if,"")


singletons_tax <- left_join(df,lineage_df, by="tax_id")


singletons_tax_genus <- 
  singletons_tax %>% 
  #bind_rows(tibble(classified= "C", seq_id="aaa", tax_id=564608, matches=as.double(84), kmer_counter="kaiju", type="kaijuresults_genbank_eukaryotes_171014_kaiju_db_85_20", rank="no", kingdom="2759", phylum="3041", class="1035538", order="13792", family="41873", genus=NA, species=NA),.) %>% 
  mutate(
    final_tax_id = ifelse(
      str_detect(rank, "species")|!is.na(species),
      paste(tax_id, kingdom, phylum, class, order, family, genus, sep="__") %>% str_remove_all(pattern = "NA__|__NA|NA") %>% str_remove("\\w+__") ,
      ifelse(is.na(genus), tax_id, genus )
    )
  ) %>% 
  select(classified, seq_id, genus, final_tax_id)

singletons_tax_genus %>% 
  filter(genus==final_tax_id) %>% 
  count() %>%
  rename("count"=n) %>% 
  bind_cols(type="genus",.) %>% 
  write_tsv(singletons_genus_stats_file)

singletons_tax_genus %>% 
  select(classified, seq_id, final_tax_id ) %>% 
  rename("tax_id"=final_tax_id) %>% 
  write_tsv(singletons_genus_file)




# singletons_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage4/taxonomy_processing/singletons_added_SGF_empty.txt"
# lineage_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage4/taxonomy_processing/singletons_added_SGF_empty_lineage.txt"
# 
# singletons_genus_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/SE_RNA/stage4/taxonomy_processing/singletons_genus.txt"
# singletons_genus_stats_file <- "/Users/pernillaericsson/Documents/medair1/apps/bio/dev_repos/parca/demo/snakemake_results_a/stats_SE_RNA/stage4/count_singletons_genus.txt"

# mutate(
  #   final_tax_id = ifelse(str_detect(rank, "species"),  
  #                         paste(tax_id, kingdom, phylum, class, order, family, genus, sep="__") %>% str_remove_all(pattern = "NA__|__NA|NA") %>% str_remove("\\w+__")  ,
  #                         ifelse(is.na(species), 
  #                               ifelse(is.na(genus), tax_id, genus ),  
  #                               paste(tax_id, kingdom, phylum, class, order, family, genus, sep="__") %>% str_remove_all(pattern = "NA__|__NA|NA") %>% str_remove("\\w+__") )
  # 
  # ) ) %>% 
  # select(seq_id, final_tax_id)

    # mutate(
  #   finaltaxid = ifelse(str_detect(rank, "species")|is.na(species),  
  #                       ifelse(is.na(genus),tax_id, genus ), 
  #                       tax_id )
  # )
  #filter(str_detect(rank, "class")) %>% pull(rank) %>% unique()
  #pull(rank) %>% unique()
  # mutate(
  #   finaltaxid = ifelse(is.na(genus)|is.na(species),tax_id, genus)
  # ) %>% 
  # filter(is.na(genus)) %>% #%>% pull(rank) %>% unique()
  # filter(str_detect(rank,"species"))
  # 
# singletons_tax %>% 
#   unite("lowest_vec", c(kingdom, phylum, class, order, family, genus), na.rm = TRUE, remove = FALSE, sep="__") %>% 
#   mutate(lowest = str_remove(lowest_vec,"\\w+__") )

#superkingdom,class,order,family,genus,species