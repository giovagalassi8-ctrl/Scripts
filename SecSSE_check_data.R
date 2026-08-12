#!/usr/bin/env Rscript

# This script validates and summarizes phylogenetic tree properties, trait coverage, and sampling fraction for a selected clade prior to SecSSE analysis. 
# It cross-references tree tip labels with the trait table to guarantee perfect matching, evaluates tree topology (e.g., ultrametric, binary),
# and calculates the true sampling fraction based on the accepted total number of extant species (obtained consulting online databases or literature).
# Depth data is the data type chosen as the default for the script, but can be changed to any other environmental data type for the species under consideration.
# It requires a Newick file, containing species names written with underscores, and the .csv file obtained after runnning the 'SecSSE_prepare_trait.R' script.

# Once you have checked the output file, if the results are corrected, you could start with the models analyses for the SecSSE approach. 


library(ape)


# Define the input treefile in .nwk format.
tree_file <- "TREEFILE.NWK"
# Definte the .csv file containing the data for the environmental parameter selected for the analysis.
# It could be used the output file of the 'SecSSE_trait_data.R' script.
trait_file <- "TRAIT_FILE.CSV"
# Define the output file name.
output_file <- "OUTPUT_DATA_CHECK_FILE.CSV"

# Define the estimated total number of accepted extant species.
# To gain this number you can use online databases (e.g., WoRMS for marine species).
total_known_species <- 610   #(Change accordingly).

# Load the phylogenetic tree.
tree <- read.tree(tree_file)
# Check if the file contains multiple trees and just in case it uses only the first one. 
if (inherits(tree, "multiPhylo")) {
  tree <- tree[[1]]
}

# Load the table containing the environmental trait data.
trait_data <- read.csv(trait_file, stringsAsFactors = FALSE, na.strings = c("NA", ""))

# Convert tree labels from underscores to spaces to match with the species name in the environment file.
# If the environment file contains species name with the underscore (e.g., Verpa_penis), skip this line.
trait_data$tree_label <- gsub(" ", "_", trait_data$species)

# Identify mismatches between the treefiles and trait table.
# Create a vector of species in the treefile but not in the trait table.
tree_not_data <- setdiff(
  tree$tip.label,
  trait_data$tree_label)
# Create a vector of species in the trait table but not in the treefile.
data_not_tree <- setdiff(
  trait_data$tree_label,
  tree$tip.label)

# Stop the process if the species lists differ.
if (length(tree_not_data) > 0 || length(data_not_tree) > 0) {
  stop("Tree and trait species do not match.")
}

# Reorder the trait table to match the tree tips.
trait_data <- trait_data[
  match(tree$tip.label, trait_data$tree_label),
]

# Extract the trait vector as numeric.
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
  na.rm = TRUE)

# Create the diagnostic summary.
check_results <- data.frame(
  number_of_tips = Ntip(tree),
  total_known_species = total_known_species,
  sampling_fraction = sampling_fraction_value,
  # Count the species classified as 'shallow' (trait = 0), explicitly ignoring missing data.
  # In this case depth is used as default parameter. If another environmental parameter is considered, change with the correct name (not shallow).
  number_shallow = sum(trait == 0, na.rm = TRUE),
  # Count the species classified as 'deep' (trait = 1), explicitly ignoring missing data.
  # In this case depth is used as default parameter. If another environmental parameter is considered, change with the correct name (not deep).
  number_deep = sum(trait == 1, na.rm = TRUE),
  # Count the total number of species in the tree that lack data.
  number_missing = sum(is.na(trait)),
  # Compute the percentage of species in the tree that possess valid, non-NA trait data.
  percentage_with_trait = mean(!is.na(trait)) * 100,
  # Record the calculated mean parameter value used as the threshold to categorize species.
  # Change with the correct environmental parameter.
  depth_threshold = depth_threshold,
  # Check if the tree is rooted.
  tree_is_rooted = is.rooted(tree),
  # Check if the tree is ultrametric.
  tree_is_ultrametric = is.ultrametric(tree),
  # Check if the tree is strictly bifurcating and contains no polytomies.
  tree_is_binary = is.binary.tree(tree),
  # Confirm that the tree object structurally contains branch length data.
  branch_lengths_present = branch_lengths_present,
  # Validate that all existing branch lengths are finite and strictly positive.
  branch_lengths_valid = branch_lengths_valid,
  # Guarantee perfect sequential matching between trait labels and tree tips.
  species_lists_match = identical(trait_data$tree_label, tree$tip.label)
)

# Save the diagnostic summary.
write.csv(check_results, output_file, row.names = FALSE)
