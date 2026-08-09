#!/usr/bin/env Rscript

# This script retrieves a time-calibrated phylogenetic tree (timetree) for any specified animal family.
# It fetches the standard topology from the Open Tree of Life (OTL) and applies a root age calibration based on user-provided 
# literature data to estimate branch lengths in millions of years (Ma).
# It requires a taxonomic clade name and the crown age of the group in Millions of years (Ma).
# This script could be used to create a time calibrated tree that could be used as 'resolve_taxonomy_marine_animals.R' input if the group of interest has no data in the timetree.org database. 


library(rotl)
library(ape)

# Set the clade of interest name (even if TARGET_FAMILY is referenced, other levels like genus can also be targeted).
TARGET_FAMILY      <- "CLADE_OF_INTEREST_NAME"  # Check if there are data for that group on OTL before.
# Change with the crown age in Ma, from the literature.
ROOT_AGE_MA        <- Ma

# Create the output file name based on the target group.
output_file <- paste0(TARGET_FAMILY, ".nwk")

# Fetching topology from Open Tree of Life.
# Match the target family name against the OTL taxonomy Reference System.
match   <- tnrs_match_names(TARGET_FAMILY, context_name = "Animals")
# Halt the execution if the family name is not recognuzed by OTL.
if (is.na(match$ott_id[1])) stop("Family not found in OTL - check spelling.")
# Extract the phylogenetic subtree associated with the identified OTL ID.
tree <- tol_subtree(ott_id = match$ott_id[1])
cat("Tips from OTL:", Ntip(tree), "\n")

# Clean uninformative tips.
# Remove tips that lack clear species-level identification to prevent issues during downstream database matches (e.g., GBIF).
# Identify tips containing "incertae sedis", "unclassified", or lacking an underscore.
to_drop <- tree$tip.label[
  grepl("incertae.sedis|unclassified", tree$tip.label, ignore.case = TRUE) |
    !grepl("_", tree$tip.label)
]
# Drop the identified uninformative tips from the tree structure.
if (length(to_drop) > 0) {
  cat("Removing", length(to_drop), "uninformative tips\n")
  tree <- drop.tip(tree, to_drop)
}
cat("Tips after cleaning:", Ntip(tree), "\n")

# Ensure the tree has enough tips left to perform a meaningful analysis.
if (Ntip(tree) < 4) stop("Fewer than 4 tips - too small for analysis.")

# Clean tip names.
# Reformat tip labels to standard binomial nomenclature for compatibility.
# Strip the OTL internal identifiers (e.g., "_ott512553") from every tip label.
tree$tip.label <- sub("_ott[0-9]+$", "", tree$tip.label)

# Some tips may still have more than two name components after OTL code removal.
# There are three possible cases:
# 1. Environmental samples / non-species labels: removed entirely (e.g. "Hormogaster_environmental_sample");
# 2. Names with "aff.", "cf.", "subsp." or years: removed entirely;
# 3. Subspecies (exactly 3 underscore-separated words): truncated to the first twp names only, so GBIF can find occurrence records under the species name. 
# If two subspecies of the same species are both in the tree, they will collapse to the same binomial.
to_remove_bad <- tree$tip.label[
  grepl("environmental_sample|_aff\\._|_cf\\._|subsp\\.|[0-9]{4}", 
        tree$tip.label, ignore.case = TRUE)
]
# Drop the invalid taxa from the tree.
if (length(to_remove_bad) > 0) {
  cat("Removing", length(to_remove_bad), "non-species tips:\n")
  print(to_remove_bad)
  tree <- drop.tip(tree, to_remove_bad)
}

# Identify trinomials (subspecies) by counting the number of underscores.
is_trinomial <- lengths(regmatches(tree$tip.label,
                                   gregexpr("_", tree$tip.label))) >= 2
# Truncate trinomials down to standard binomials.
if (any(is_trinomial)) {
  cat("\nTruncating", sum(is_trinomial), "subspecific names to binomial:\n")
  print(tree$tip.label[is_trinomial])
  # Capture the first two underscore-separated words, ignore the rest.
  tree$tip.label[is_trinomial] <- sub("^([^_]+_[^_]+)_.*$", "\\1",
                                      tree$tip.label[is_trinomial])
}

# Check for duplicate binomials (two subspecies collapsed to the same name).
dupes <- tree$tip.label[duplicated(tree$tip.label)]
if (length(dupes) > 0) {
  cat("\nWARNING: duplicate tip names after truncation (keeping first occurrence):\n")
  print(unique(dupes))
  # Identify the indices of the duplicate tips and drop them, keeping only the first.
  to_drop_dupes <- which(duplicated(tree$tip.label))
  tree <- drop.tip(tree, tree$tip.label[to_drop_dupes])
}
cat("Final tip count:", Ntip(tree), "\n")

# Time calibration process: convert the strict topology into a chronogram (time-scaled tree).
# Resolve any polytomies (nodes with >2 descendants) randomly into dichotomies.
if (!is.binary(tree)) tree <- multi2di(tree)
# The chronos() function requires a fully bifurcating (binary) tree.
# Assign a uniform arbitrary branch length of 1 to all branches.
# This is necessary because OTL trees lack branch lengths, and chronos() needs initial non-zero numeric values to begin optimization.
tree$edge.length <- rep(1, nrow(tree$edge))
# Create a calibration object anchoring the root node (Ntip + 1).
# The calibration bounds are set to ±10% of the provided ROOT_AGE_MA.
calib       <- makeChronosCalib(tree, node = Ntip(tree) + 1,
                                age.min = ROOT_AGE_MA * 0.9,
                                age.max = ROOT_AGE_MA * 1.1)
# Apply the penalized likelihood method to estimate branch lengths in time.
timetree    <- chronos(tree, calibration = calib, lambda = 1)

# Verify that the resulting tree is ultrametric (all tips end at time 0).
if (!is.ultrametric(timetree)) stop("Tree is not ultrametric - chronos() failed.")


# Export the time-calibrated tree to a standard Newick file format
write.tree(timetree, file = output_file)
cat("Tips:", Ntip(timetree), "| Root age:", 
    round(max(branching.times(timetree)), 1), "Ma\n")

# Generate the visual representation of the timetree to check if it is correct.
# Plot the phylogenetic tree.
plot(timetree,
     #Tip size scales dynamically based on tree size.
     cex      = min(0.8, 40 / Ntip(timetree)),
     # Set the plot title based on the name of the family.
     main     = TARGET_FAMILY,
     no.margin = FALSE)
# Add a time axis (in Ma) at the bottom of the plot.
axisPhylo()
# Add a descriptive subtitle displaying the root age.
mtext(paste0("Root age: ", ROOT_AGE_MA, " Ma"),
      side = 1, 
      line = 3, 
      cex = 0.8, 
      col = "grey30")
