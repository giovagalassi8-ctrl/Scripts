
# Purpose:
# Summarise tree properties, trait coverage and sampling fraction.
#
# Required trait CSV columns:
#   species
#   Mean_Depth
#   Depth_magnitude
#   trait
#
# The trait table must contain the same species as the phylogenetic tree.
#
# Output:
#   analyses/SSE/Terebridae_data_check.csv
# =============================================================================


library(ape)


# Define input and output files.
tree_file <- "data/Mollusca/Terebridae/Terebridae.nwk"
trait_file <- "Terebridae_depth_trait.csv"
output_file <- "Terebridae_data_check.csv"


# Define the estimated total number of accepted extant species.
total_known_species <- 610


# Load the phylogenetic tree.
tree <- read.tree(tree_file)

if (inherits(tree, "multiPhylo")) {
  tree <- tree[[1]]
}


# Load the trait table.
trait_data <- read.csv(
  trait_file,
  stringsAsFactors = FALSE,
  na.strings = c("NA", "")
)


# Convert species names to the tree-label format.
trait_data$tree_label <- gsub(
  " ",
  "_",
  trait_data$species
)


# Identify mismatches between the tree and trait table.
tree_not_data <- setdiff(
  tree$tip.label,
  trait_data$tree_label
)

data_not_tree <- setdiff(
  trait_data$tree_label,
  tree$tip.label
)


# Stop if the species lists differ.
if (length(tree_not_data) > 0 || length(data_not_tree) > 0) {
  stop("Tree and trait species do not match.")
}


# Reorder the trait table to match the tree tips.
trait_data <- trait_data[
  match(tree$tip.label, trait_data$tree_label),
]


# Extract the trait vector.
trait <- as.numeric(trait_data$trait)


# Calculate the sampling fraction.
sampling_fraction_value <- Ntip(tree) / total_known_species


# Check whether branch lengths are available and valid.
branch_lengths_present <- !is.null(tree$edge.length)

branch_lengths_valid <- (
  branch_lengths_present &&
    all(is.finite(tree$edge.length)) &&
    all(tree$edge.length > 0)
)


# Calculate the depth threshold used for classification.
depth_threshold <- mean(
  trait_data$Depth_magnitude,
  na.rm = TRUE
)


# Create the diagnostic summary.
check_results <- data.frame(
  number_of_tips = Ntip(tree),
  total_known_species = total_known_species,
  sampling_fraction = sampling_fraction_value,
  number_shallow = sum(trait == 0, na.rm = TRUE),
  number_deep = sum(trait == 1, na.rm = TRUE),
  number_missing = sum(is.na(trait)),
  percentage_with_trait = mean(!is.na(trait)) * 100,
  depth_threshold = depth_threshold,
  tree_is_rooted = is.rooted(tree),
  tree_is_ultrametric = is.ultrametric(tree),
  tree_is_binary = is.binary.tree(tree),
  branch_lengths_present = branch_lengths_present,
  branch_lengths_valid = branch_lengths_valid,
  species_lists_match = identical(
    trait_data$tree_label,
    tree$tip.label
  )
)


# Save the diagnostic summary.
write.csv(
  check_results,
  output_file,
  row.names = FALSE
)
