# !/usr/bin/env Rscript

# This script takes a phylogenetic tree (in Newick format) and looks up each species name in WoRMS and GBIF databases to find the currently accepted scientific name.
# Sometimes the tip labels often use older names or synonyms that no longer match the name used in the database when the datas were added.
# This could be useful if you have to download geographic occurrence data for the species in a phylogenetic tree.
# The final result of this script could be used to run the "geographic_coordinates_gbif.R" .

# IMPORTANT: Both WoRMS and GBIF use approximate name matching. 
# This means that for a name with no exact match in the database, the script will still return its closest guess rather then falling outright (the guess is not guaranteed to be correct).
# A manual read-through of the final CSV is still strongly recommended, because the automated name matching can occasionally pick the wrong species.


library(ape)
library(worrms)
library(rgbif)
library(stringr)

# Import the phylogenetic tree (mainly a time tree, downloaded from timetree.org) in Newick format (change accordingly).
tree <- read.tree("NEWICK_TREE")
# Create an object that includes only the tips labels.
tip_names <- tree$tip.label
# Replace the underscore in the species name labels with a space (e.g., Octopus_vulgaris -> Octopus vulgaris).
# Taxonomic database usually expect a binominal species name with a space.
tip_species <- str_replace_all(tip_names, "_", " ")

# Start the research for the species in the WoRMS database.
# Set as input a single species name at a time.
resolve_worms <- function(sp_name) {
  res <- tryCatch(
    # Performs a fuzzy matching against the WoRMS database and is tolerant of minor spelling differences. 
    # It also resolves synonyms automatically, returning the currently accepted name (valid_name) together with its unique WoRMS identifier (AphiaID).
    wm_records_taxamatch(
      name = sp_name,
      # In this case, the research is made only for marine animals, so if there is a terrestrial synonym it will not be considered.
      marine_only = TRUE
    ),
    error = function(e) NULL
  )
  # If WoRMS returns nothing (no match at all), flag that name as NOT_FOUND so it can be picked up by the GBIF fallback step below.
  if (is.null(res) || length(res) == 0 || nrow(res[[1]]) == 0) {
    return(data.frame(
      queried_name = sp_name,
      status = "NOT_FOUND",
      valid_name = NA, 
      valid_AphiaID = NA,
      n_candidate_matches = 0
    ))
  }
  
  # Record how many candidates row WoRMS returned in total before discarding all but the first.
  n_candidates <- nrow(res[[1]])
  # Take the first (best) match returned by the fuzzy search. The first row in the results of the research is supposed to be the best, but it is better to check it.
  rec <- res[[1]][1, ]
  # Build a one-row summary for this species, keeping only the fields we need.
  data.frame(
    queried_name = sp_name, # Write the name that appears on the tree.
    status = rec$status, # Write the WoRMS status of the species (e.g., accepted, synonym, ...).
    valid_name = rec$valid_name, # Write the currently accepted name. 
    valid_AphiaID = rec$valid_AphiaID, # Write the WoRMS taxon identifier.
    n_candidate_matches = n_candidates # Write the number of possible WoRMS matches for that name (if >1, only the top ranked one was kept, and this row deserves an extra check).
  )
}

# Apply the WoRMS lookup to every tip in the tree.
worms_results <- do.call(rbind, lapply(tip_species, resolve_worms))
# Flag matches that are NOT an exact string match between the queried name and the returned valid_name.
# It narrows down which rows deserve a closer look during manual review.
worms_results$exact_match <- mapply(function(queried, valid) {
  if (is.na(valid)) return(NA)
  identical(queried, valid)
}, worms_results$queried_name, worms_results$valid_name)
 
cat("\nRows where the resolved name differs from the queried name:\n")
print(worms_results[!isTRUE(worms_results$exact_match) & !is.na(worms_results$valid_name), ])

cat("\nRows where WoRMS found more than one possible match:\n")
print(worms_results[!is.na(worms_results$n_candidate_matches) & worms_results$n_candidate_matches > 1, ])

# Print the complete table.
print(worms_results)

# For any species names WoRMS could not resolve, the species will be searched on GBIF instead.
# GBIF sometimes recognizes names that WoRMS does not, in particular recently described species.
# Create the list of species that WoRMS could not match at all.
not_found <- worms_results$queried_name[worms_results$status == "NOT_FOUND"]

# Try the GBIF research only if there is at least one unresolved name.
# Skip this whole block when WoRMS has already resolved every tip.
if (length(not_found) > 0) {
  cat("\nSpecies not found in WoRMS, attempting GBIF fallback:\n")
  print(not_found)
  # For each unresolved name, the research will be done using the GBIF database (in this case it is not limited to marine animals).
  gbif_fallback <- lapply(not_found, function(nm) {
    res <- name_backbone(name = nm, rank = "species")
    # Build a one-row summary for this species, keeping only the fields we need.
    data.frame(
      # Each ifelse() substitutes NA if GBIF returned nothing for that field, so every row has the same structure even with no match.
      queried_name = nm, # Write the name as it appeared in the tree.
      gbif_matchType = ifelse(is.null(res$matchType), NA, res$matchType), # Write the level of the GBIF match (e.g., exact, fuzzy, ...).
      gbif_status = ifelse(is.null(res$status), NA, res$status), # Write the GBIF taxonomic status (e.g., accepted, synonym, ...).
      gbif_accepted_name = ifelse(is.null(res$species), NA, res$species), # Write the name GBIF currently considers valid for this species.
      gbif_usageKey = ifelse(is.null(res$usageKey), NA, res$usageKey) # Write the numeric ID of GBIF for this taxon.
    )
  })
  # Stack all the one-rows results into a single dataframe.
  gbif_fallback <- do.call(rbind, gbif_fallback)
  # Show the results on screen.
  print(gbif_fallback)

  # Save as separate file from the WoRMS results, since the column names are different (change the file name accordingly).
  write.csv(gbif_fallback, "GBIF_TAXONOMY.CSV", row.names = FALSE)
}

# Save the final reference file (change the file name accordingly).
write.csv(worms_results, "WORMS_TAXONOMY.CSV", row.names = FALSE)

cat("Resolved via WoRMS:", sum(worms_results$status != "NOT_FOUND"), "/", length(tip_species), "\n")
cat("Still unresolved (NOT_FOUND):", sum(worms_results$status == "NOT_FOUND"), "\n")
