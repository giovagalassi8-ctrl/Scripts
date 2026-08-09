#!/usr/bin/env Rscript

# This script prepare the binary trait data used for the SecSSE macroevolutionary analysis.
# This script extracts species from a phylogenetic tree and matches them with environmental depth data.
# Depth data is the data type chosen as the default for the script, but can be changed to any other environmental data type for the species under consideration.
# A classification threshold is dynamically calculated based on the mean depth of the dataset, categorizing the species into a binary state:
#   - 0: Shallow-water species (Depth <= Mean Depth)
#   - 1: Deep-water species (Depth > Mean Depth)
#   - NA: Missing depth information
# It requires a Newick file, containing species names written with underscores, and a .csv file containing three specific columns: Clade, species and Mean_Depth

# The output of this script could be used as input for the downstream SecSSE script 'SecSSE_check_data.R'.


library(ape)
library(dplyr)


# Define input and output files.
tree_file <- "data/Mollusca/Terebridae/Terebridae.nwk"
environment_file <- "analyses/02_Tip_rates/marine_environment.csv"
output_file <- "Terebridae_depth_trait.csv"


# Load the phylogenetic tree.
tree <- read.tree(tree_file)

if (inherits(tree, "multiPhylo")) {
  tree <- tree[[1]]
}


# Convert tree labels from underscores to spaces.
species_names <- gsub("_", " ", tree$tip.label)


# Load the environmental data.
environment_data <- read.csv(
  environment_file,
  stringsAsFactors = FALSE
)


# Keep only the focal clade and the required variables.
environment_data <- environment_data %>%
  filter(Clade == "Terebridae") %>%
  select(species, Mean_Depth)


# Convert depth values to numeric.
environment_data$Mean_Depth <- as.numeric(
  environment_data$Mean_Depth
)


# Create the complete species list from the tree.
trait_data <- data.frame(
  species = species_names,
  stringsAsFactors = FALSE
)


# Add depth information to the tree species.
trait_data <- trait_data %>%
  left_join(
    environment_data,
    by = "species"
  )


# Convert depth values to absolute magnitude.
trait_data$Depth_magnitude <- abs(
  trait_data$Mean_Depth
)


# Calculate the mean depth used as classification threshold.
depth_threshold <- mean(
  trait_data$Depth_magnitude,
  na.rm = TRUE
)


# Classify species as shallow or deep.
trait_data <- trait_data %>%
  mutate(
    trait = case_when(
      is.na(Depth_magnitude) ~ NA_real_,
      Depth_magnitude <= depth_threshold ~ 0,
      Depth_magnitude > depth_threshold ~ 1
    )
  )


# Save the trait table.
write.csv(
  trait_data,
  output_file,
  row.names = FALSE,
  na = "NA"
)
