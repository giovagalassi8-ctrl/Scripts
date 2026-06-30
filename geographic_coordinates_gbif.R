# !/usr/bin/env Rscript

# This script downloads geographic occurrence records from the GBIF database, for a specific list of marine species.
# It requires teh WORMS_TAXONOMY file obtained by the 'resolve_taxonomy_marine.R', from which it gather, 
# filter and standardize geographic occurrence records from the GBIF database.
# The output could be used to run the 'geographic_overlap.R' script.

library(rgbif)
library(dplyr)

# Import the csv file obtained by the "resolve_taxonomy_marine.R" containing the currently accepted scientific names (change accordingly).
taxonomy <- read.csv("WORMS_TAXONOMY.CSV", stringsAsFactors = FALSE)

# Build the actual name to send to GBIF for each species.
# In this case, the resolved 'valid_name' is preferred if one was found; 
# if a species could not be resolved at all, the original tree tip name are selected.
taxonomy$name_to_query <- ifelse(!is.na(taxonomy$valid_name) & taxonomy$valid_name != "",
                               taxonomy$valid_name, taxonomy$queried_name)

# Shows the final list of names about to be queried, for a quick check.
print(lookup$name_to_query)

# This function takes a single species name and returns its GBIF occurrence data as a data frame (or NULL).
download_gbif <- function(species_name) {
  res <- tryCatch(
    # queries GBIF for occurrence records matching this scientific name.
    occ_search(
      scientificName = species_name,
      # Discards any record without usable latitude and longitude coordinates.
      hasCoordinate = TRUE,
      # Limit the number of records to 2000 per species to prevent excessive download times (change desired).
      limit = 2000),
    # Define what happens if an error occurs during the occ_search().
    error = function(e) {
      cat("  Error for", species_name, ":", conditionMessage(e), "\n")
      NULL
    }
  )
  
  # Check if the result is NULL, if the data field is missing, or if there are zero rows returned.
  if (is.null(res) || is.null(res$data) || nrow(res$data) == 0) {
    cat("  No occurrences found for", species_name, "\n")
    return(NULL)
  }
  # Return the actual occurrence dataframe.
  res$data
}

# Apply the download_gbif function to every species name in the name_to_query column.
all_occurrences <- lapply(taxonomy$name_to_query, download_gbif)
# Filter out any NULL elements from the list (species that had no data or caused an error)
all_occurrences <- all_occurrences[!sapply(all_occurrences, is.null)]


# Define the essential columns to keep in order to avoid type conflicts across different species dataframes (change if necessary).
essential_cols <- c("scientificName", "species", "decimalLatitude",
                    "decimalLongitude", "depth", "country", "year")

# Iterate over the list of dataframes to extract only the essential columns.
combined <- lapply(all_occurrences, function(df) {
  # Identify which of the essential columns are actually present in the current dataframe.
  cols_present <- intersect(essential_cols, colnames(df))
  # Subset the dataframe to keep only those present columns, preventing reduction to a vector.
  df[, cols_present, drop = FALSE]
})
# Bind all the individual dataframes together row by row into one master dataframe.
combined_df <- bind_rows(combined)

# Print the total number of occurrence records successfully downloaded.
cat("\nTotal occurrences downloaded:", nrow(combined_df), "\n")
# Print the number of species that have at least one occurrence against the total number of queried species.
cat("Species with at least one occurrence:", length(unique(combined_df$species)), "/", nrow(taxonomy), "\n")

# Count the number of occurrences for each species and sort them in ascending order.
occ_per_species <- combined_df %>% count(species, name = "n_occurrences") %>% arrange(n_occurrences)
print(occ_per_species)

# Define the minimum threshold of occurrences considered reliable for downstream analysis.
MIN_OCC <- 5
# Identify species that fall below this minimum occurrence threshold.
sp_pochi_dati <- occ_per_species$species[occ_per_species$n_occurrences < MIN_OCC]
# Check if there is at least one species with too few occurrences.
if (length(sp_pochi_dati) > 0) {
  # Warn the user that these species have very few records, making geographical overlap unreliable.
  cat("\nWARNING - species with fewer than", MIN_OCC, "occurrences (overlap unreliable):\n")
  print(sp_pochi_dati)
}

# Identify species that were queried but are completely missing from the final GBIF dataframe.
missing_species <- setdiff(taxonomy$name_to_query, unique(combined_df$species))
# Check if there are any completely missing species.
if (length(missing_species) > 0) {
  cat("\nWARNING - species with ZERO GBIF occurrences:\n")
  print(missing_species)
}

# Write the combined dataframe to a CSV file without including row names (change file name accordingly).
# This file could be used to run the 'geographic_overlap.R' script.
write.csv(combined_df, "COORDINATES.CSV", row.names = FALSE)
