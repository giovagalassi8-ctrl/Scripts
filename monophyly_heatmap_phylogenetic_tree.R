#!/usr/bin/env Rscript

# This script generates a phylogenetic tree annotated with a heatmap.
# It visualizes the monophyly status (true/false) of different taxonomic groups across multiple alternative trees.

# REQUIRED FILES:
# 1. A phylogenetic tree file in Newick format;
# 2. A CSV file containing the monophyly results from a gotree analysis.

# Load necessary libraries
library(ape)
library(ggtree)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)

# Load the phylogenetic tree from a Newick file.
tree <- read.tree("03_Gotree/03_treefile/Spiralia_semplified_tree.nwk")

# Load the monophyly status data from a CSV file.
monophyly <- read.csv("03_Gotree/00_stats/monophyly_results_gotree.csv", stringsAsFactors = FALSE)

# The two following steps are necessary to have a file that can be used in creating the graph with the next function.
# The ggplot2 library strictly requires that data be organized in a "long" format.
monophyly_long <- pivot_longer(monophyly, cols = -Tree, names_to = "Group", values_to = "Status")
# The ggtree library needs another conversion, this time into wide format, using "Tree" names as columns and "Group" as rows.
monophyly_wide <- pivot_wider(monophyly_long, names_from = Tree, values_from = Status)

# Convert the "Group" column into row names, which is required by the gheatmap function to map to the tree.
heatmap_data <- column_to_rownames(monophyly_wide, var = "Group")

# Build the phylogenetic tree plot.
phylo_tree <- ggtree(tree) +
  # Add tip labels and align them with dotted lines
  geom_tiplab(align = TRUE, linesize = 0.5, offset = 0.1, size = 4) +
  # A clean theme for trees
  theme_tree() 

# Add the heatmap to the phylogenetic tree
heatmap_plot <- gheatmap(phylo_tree, heatmap_data,
                         offset = 4,           # Distance between tip labels and the first heatmap column.
                         width = 4.6,            # Total width of the heatmap relative to the tree.
                         colnames = FALSE,        # Don't show the names of the trees with FALSE (digit TRUE if you want to see them).
                         color = "black",        # Draws black borders around tiles, simulating separated squares.
                         font.size = 3) +
 # Define custom fill colors for the monophyly status.
  scale_fill_manual(values = c("false" = "white", "true" = "darkgreen", 
                               "FALSE" = "white", "TRUE" = "darkgreen"),
                    # Ensures missing data (taxa without a match between the tree and the table) does not get a background color.
                    na.value = "transparent",
                    name = "Monophyly Status") +
  # Customize the legend position and the title appearance.
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "MS90")


# Fine-tune the heatmap aesthetics by directly modifying the ggplot layers
# Find the index of the heatmap layer (it is typically the last layer added to the plot)
heatmap_layer <- length(heatmap_plot$layers)

# Extract the data directly from the heatmap layer to manipulate it
layer_data <- heatmap_plot$layers[[heatmap_layer]]$data

# Remove NA values completely
# This prevents the plot from drawing empty black-bordered boxes for missing data
heatmap_plot$layers[[heatmap_layer]]$data <- subset(layer_data, !is.na(value))

# Calculate the exact x-axis spacing between columns
# Find the unique x-coordinates for the heatmap columns and calculate the distance between the first two
x_coords <- sort(unique(layer_data$x))
col_spacing <- x_coords[2] - x_coords[1]

# Apply dynamic shrinking for both rows and columns to create separated square tiles
# Shrink the height of the tiles (the y-axis spacing is always exactly 1 integer unit in ggtree)
heatmap_plot$layers[[heatmap_layer]]$aes_params$height <- 0.8
# Shrink the width of the tiles to 80% of the calculated dynamic column spacing
heatmap_plot$layers[[heatmap_layer]]$aes_params$width <- col_spacing * 0.8


# Save the final plot to a PDF file
ggsave(file = "phylotree_monophyly_heatmap.pdf", plot = heatmap_plot, width = 10, height = 7, dpi = 300)
