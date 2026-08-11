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


# Define the input treefile in .nwk format (change accordingly).
tree_file <- "TREEFILE.NWK"
# Define the .csv file containing the environmental data (the script has depth as default parameter analysed).
# Change accordingly (You can use the marine_environment.csv or the terrestrial_environment.csv obtained using 'extraction_marine_environment.R' or 'extraction_terrestrial_environment.R', respectively).
environment_file <- "ENVIRONMENT.CSV"
# Define the output file name.
output_file <- "OUTPUT_TRAIT.CSV"


# Load the phylogenetic tree.
tree <- read.tree(tree_file)
# Check if the file contains multiple trees and just in case it uses only the first one. 
if (inherits(tree, "multiPhylo")) {
  tree <- tree[[1]]
}

# Convert tree labels from underscores to spaces to match with the species name in the environment file.
# If the environment file contains species name with the underscore (e.g., Verpa_penis), skip this line.
species_names <- gsub("_", " ", tree$tip.label)

# Load the environmental data.
environment_data <- read.csv(environment_file, stringsAsFactors = FALSE)

# Keep only the focal clade and the required variables.
environment_data <- environment_data %>%
  # Change with the correct clade name.
  filter(Clade == "CLADE_NAME") %>%
  # Change with the name of the file's column containing the parameter of interest. Depth is the default parameter, but you can use any other environmental condition.
  # If you change the following line, you have to change also the downstream command line with the correct name.
  select(species, Mean_Depth)

# Convert depth values to numeric.
# If you use a different parameter, change the name accordingly with the previous line.
environment_data$Mean_Depth <- as.numeric(environment_data$Mean_Depth)

# Create the complete species list from the tree.
trait_data <- data.frame(
  species = species_names,
  stringsAsFactors = FALSE)

# Add depth information to the tree species.
# Change with the correct environmetal parameter used.
trait_data <- trait_data %>%
  left_join(
    environment_data,
    by = "species")

# Convert depth values to absolute magnitude.
# Change with the correct environmetal parameter used.
trait_data$Depth_magnitude <- abs(trait_data$Mean_Depth)

# Calculate the mean depth used as classification threshold.
depth_threshold <- mean(
  trait_data$Depth_magnitude,
  na.rm = TRUE)

# Classify species as shallow or deep based on their values relative to the mean depth previously calculated.
# Change this command line if a different parameter was used (also the structure if necessary).
trait_data <- trait_data %>%
  mutate(
    trait = case_when(
      is.na(Depth_magnitude) ~ NA_real_,
      Depth_magnitude <= depth_threshold ~ 0,
      Depth_magnitude > depth_threshold ~ 1))

# Save the trait table.
write.csv(trait_data, output_file, row.names = FALSE, na = "NA")
