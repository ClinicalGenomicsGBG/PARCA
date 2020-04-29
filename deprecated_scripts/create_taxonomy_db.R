#
#
suppressPackageStartupMessages({ 
  library(taxonomizr) 
  } )


names_dmp_file <- snakemake@input[['names_dmp']]
nodes_dmp_file <- snakemake@input[['nodes_dmp']]
sql_db_file <- snakemake@output[['sql_db']]

read.names.sql(names_dmp_file, sql_db_file)

read.nodes.sql(nodes_dmp_file,sql_db_file)