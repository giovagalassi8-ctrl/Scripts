# !/usr/bin/env/ Rscript

# This script merges a BUSCO summary table with a species-to-clade mapping file.
# It allows you to evaluate the quality of BUSCO (Benchmarking Universal Single-Copy Orthologs) genomes, 
# which measures the completeness of a genomic assembly or transcriptome.
# The script organizes the results of multiple species into a single structured table. 
# By adding clade information to each species, it becomes easy to compare the quality of assemblies between different taxonomic groups


library(dplyr)
library(readr)

# Import the .tsv file (change with the correct name).
busco <- read_tsv("01_all_short_summaries_table.tsv")

# Import the .tsv file containing the species in the first column and the corresponding clade in the second one. 
clades <- read_tsv("species_to_clade.tsv")

# Merge the two files.
# 'left_join()' merges two or multiple datasets sharing common variables and keeping all rows from the first dataset intact.
merged <- busco %>%
  left_join(clades, by = "sample")

# OPTIONAL: check if there are species without a clade.
missing <- merged %>% filter(is.na(clade))
if(nrow(missing) > 0) {
  cat("ATTENTION: Species without an assigned clade:\n")
  print(missing$sample)
}

# sort first by clade, then by the name of the species.
merged_sorted <- merged %>%
  arrange(clade, sample)

# Export the new file.
# Change the name of the output file as you want.
write_tsv(merged_sorted, "01_all_short_summaries_with_clade.tsv")

cat("File created: 01_all_short_summaries_with_clade.tsv\n")
