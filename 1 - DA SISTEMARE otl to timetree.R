# =============================================================================
# OTL to Timetree
# Gets a time-calibrated phylogenetic tree for any animal family using:
#   - Open Tree of Life for the topology
#   - A root age from the literature for the time calibration
#
# INSTRUCTIONS:
#   1. Change TARGET_FAMILY, ROOT_AGE_MA and CALIBRATION_SOURCE below
#   2. Run the script
#   3. Use the output .nwk file in the geographic overlap and LTT/gamma scripts
# =============================================================================

library(rotl)
library(ape)

# --- SETTINGS (change these for each new family) ----------------------------

TARGET_FAMILY      <- "Hormogastridae"
ROOT_AGE_MA        <- 82        # crown age in Ma, from the literature
CALIBRATION_SOURCE <- "Novo et al. 2011 (67-97 Ma)"

# ----------------------------------------------------------------------------

output_file <- paste0(TARGET_FAMILY, "_timetree.nwk")

# --- GET TOPOLOGY FROM OTL --------------------------------------------------

cat("Fetching topology from Open Tree of Life...\n")
match   <- tnrs_match_names(TARGET_FAMILY, context_name = "Animals")
if (is.na(match$ott_id[1])) stop("Family not found in OTL - check spelling.")
tree    <- tol_subtree(ott_id = match$ott_id[1])
cat("Tips from OTL:", Ntip(tree), "\n")

# --- CLEAN UNINFORMATIVE TIPS -----------------------------------------------
# Removes: incertae sedis / unclassified labels, and genus-only tips
# (no underscore = no species epithet, causes problems in GBIF downloads)

to_drop <- tree$tip.label[
  grepl("incertae.sedis|unclassified", tree$tip.label, ignore.case = TRUE) |
    !grepl("_", tree$tip.label)
]
if (length(to_drop) > 0) {
  cat("Removing", length(to_drop), "uninformative tips\n")
  tree <- drop.tip(tree, to_drop)
}
cat("Tips after cleaning:", Ntip(tree), "\n")
if (Ntip(tree) < 4) stop("Fewer than 4 tips - too small for analysis.")

# --- CLEAN TIP NAMES --------------------------------------------------------
# OTL appends an internal identifier to every tip label (e.g. "_ott512553").
# This strips everything from "_ott" onward, leaving only the binomial name.

tree$tip.label <- sub("_ott[0-9]+$", "", tree$tip.label)
# ^ strips the OTL internal identifier (e.g. "_ott512553") from every tip,
#   leaving only the species name

# Some tips may still have more than two name components after OTL code
# removal. We handle three cases:
#
#   1. Environmental samples / non-species labels: removed entirely
#      (e.g. "Hormogaster_environmental_sample" - not a real taxon)
#   2. Names with "aff.", "cf.", "subsp." or years: removed entirely
#      (these are not formally described species and won't match any database)
#   3. Subspecies (exactly 3 underscore-separated words): truncated to the
#      binomial (first two words only), so GBIF can find occurrence records
#      under the species name; note that if two subspecies of the same species
#      are both in the tree, they will collapse to the same binomial - you
#      will need to keep only one of them (the script flags this below)

to_remove_bad <- tree$tip.label[
  grepl("environmental_sample|_aff\\._|_cf\\._|subsp\\.|[0-9]{4}", 
        tree$tip.label, ignore.case = TRUE)
]
if (length(to_remove_bad) > 0) {
  cat("Removing", length(to_remove_bad), "non-species tips:\n")
  print(to_remove_bad)
  tree <- drop.tip(tree, to_remove_bad)
}

# Truncate subspecific names to binomial
is_trinomial <- lengths(regmatches(tree$tip.label,
                                   gregexpr("_", tree$tip.label))) >= 2
if (any(is_trinomial)) {
  cat("\nTruncating", sum(is_trinomial), "subspecific names to binomial:\n")
  print(tree$tip.label[is_trinomial])
  tree$tip.label[is_trinomial] <- sub("^([^_]+_[^_]+)_.*$", "\\1",
                                      tree$tip.label[is_trinomial])
}

# Check for duplicate binomials (two subspecies collapsed to the same name)
dupes <- tree$tip.label[duplicated(tree$tip.label)]
if (length(dupes) > 0) {
  cat("\nWARNING: duplicate tip names after truncation (keeping first occurrence):\n")
  print(unique(dupes))
  # keep only the first occurrence of each duplicated name
  to_drop_dupes <- which(duplicated(tree$tip.label))
  tree <- drop.tip(tree, tree$tip.label[to_drop_dupes])
}

cat("Final tip count:", Ntip(tree), "\n")

# --- TIME CALIBRATION -------------------------------------------------------

if (!is.binary(tree)) tree <- multi2di(tree)
# ^ chronos() requires a fully binary tree; multi2di() resolves polytomies

tree$edge.length <- rep(1, nrow(tree$edge))
# ^ chronos() requires branch lengths as a starting point for optimisation;
#   OTL trees have topology only (no branch lengths), so we assign a
#   uniform value of 1 to every branch - chronos() will overwrite these
#   with time-scaled values constrained by the calibration above

calib       <- makeChronosCalib(tree, node = Ntip(tree) + 1,
                                age.min = ROOT_AGE_MA * 0.9,
                                age.max = ROOT_AGE_MA * 1.1)
timetree    <- chronos(tree, calibration = calib, lambda = 1)

if (!is.ultrametric(timetree)) stop("Tree is not ultrametric - chronos() failed.")

# --- SAVE -------------------------------------------------------------------

write.tree(timetree, file = output_file)
cat("Saved:", output_file, "\n")
cat("Tips:", Ntip(timetree), "| Root age:", 
    round(max(branching.times(timetree)), 1), "Ma\n")

# --- PLOT -------------------------------------------------------------------

plot(timetree,
     cex      = min(0.8, 40 / Ntip(timetree)),
     main     = TARGET_FAMILY,
     no.margin = FALSE)
axisPhylo()
mtext(paste0("Root age: ", ROOT_AGE_MA, " Ma  |  Source: ", CALIBRATION_SOURCE),
      side = 1, line = 3, cex = 0.8, col = "grey30")
# ^ adds a subtitle at the bottom of the plot with the calibration info
