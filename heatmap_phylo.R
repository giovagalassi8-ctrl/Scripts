library(ape)
library(ggtree)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)

tree <- read.tree("03_Gotree/03_treefile/Spiralia_semplified_tree.nwk")

monophyly <- read.csv("03_Gotree/00_stats/monophyly_results_gotree.csv", stringsAsFactors = FALSE)

monophyly_long <- pivot_longer(monophyly, cols = -Tree, names_to = "Group", values_to = "Status")

monophyly_wide <- pivot_wider(monophyly_long, names_from = Tree, values_from = Status)

heatmap_data <- column_to_rownames(monophyly_wide, var = "Group")

phylo_tree <- ggtree(tree) +
  # Add tip labels and align them with dotted lines
  geom_tiplab(align = TRUE, linesize = 0.5, offset = 0.1, size = 4) +
  # A clean theme for trees
  theme_tree() 

heatmap_plot <- gheatmap(phylo_tree, heatmap_data,
                         offset = 4,           # Distance between tip labels and the first heatmap column
                         width = 4.6,            # Total width of the heatmap relative to the tree
                         colnames = FALSE,        # Show the names of the trees (from the original 'Tree' column)
                         color = "black",        # Draws white borders around tiles, simulating separated squares
                         font.size = 3) +
  # Set the fill colors and handle missing data
  scale_fill_manual(values = c("false" = "white", "true" = "darkgreen", 
                               "FALSE" = "white", "TRUE" = "darkgreen"),
                    na.value = "transparent",    # CRITICAL: Makes missing taxa completely blank (no squares)
                    name = "Monophyly Status") +
  # Final aesthetic touches
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "MS90")

# TRICK: Remove NA values directly from the heatmap layer data 
# This prevents ggplot from drawing a black border around transparent/missing squares
heatmap_layer <- length(heatmap_plot$layers)
heatmap_plot$layers[[heatmap_layer]]$data <- subset(heatmap_plot$layers[[heatmap_layer]]$data, !is.na(value))

# 6. Save the final plot to a PDF file
ggsave(file = "phylotree_monophyly_heatmap.pdf", plot = heatmap_plot, width = 10, height = 7, dpi = 300)