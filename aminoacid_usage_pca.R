#!/usr/bin/env/ Rscript

# This script performs a Preincipal Component Analysis using the amino acid frequencies (previously normalized) of the species under study,
# allowing you to visualize and identify grouping patterns.
# It is important to associate the PCA with a statistical test (in some papers about aminoacidic frequencies was used a PERMANOVA).

# USAGE:
# [Rstudio] source(aminoacid_usage_pca.R)

library(ggplot2)
library(dplyr)
library(ggrepel)

# Import dataframe containing amino acid frequencies.
data  <- read.csv("MS80_gappyout_rcv_filtered-normalized.csv", header = T)

# Define the 20 standard amino acids (single-letter codes) to use as PCA input.
aa  <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

# Run PCA on the amino acid frequency matrix.
pca <- prcomp(data[, aa], 
              center = TRUE,
              scale. = TRUE)

# Extract the percentage of total variance explained by each principal component.
# importance[2, ] is the "Proportion of Variance" row from the PCA summary.
var <- round(summary(pca)$importance[2, ] * 100, 2)

# Extract the PCA scores (coordinates of each observation in PC space),
# and attach the grouping variable so points can be colored by group.
scores        <- as.data.frame(pca$x)
scores$Group  <- data$Alignment_name

# Compute the centroid (mean PC1 and PC2) for each group.
# Centroids will be plotted as filled squares to mark the group center.
centroids <- scores %>% 
  group_by(Group) %>% 
  summarise(PC1 = mean(PC1), PC2 = mean(PC2))

# Compute the convex hull for each group, i.e. the smallest convex polygon
# that encloses all points of that group.
hulls <- scores %>% 
  group_by(Group) %>%
  # chull() returns the row indices of the boundary points in counter-clockwise order.
  # slice() selects them.
  slice(chull(PC1, PC2))

# Define a color vector for the selected groups.
group_colors <- c(
  "Macrodasyida"   = "#19647EFF",
  "Chaetonotida"   = "#BED4E9FF",
  "Platyhelmintes" = "#FDB927FF"
)

# Creates the plot.
ggplot(scores, aes(x = PC1,
                   y = PC2, 
                   color = Group,
                   fill = Group)) +
  # Draw the convex hull polygon for each group.
  geom_polygon(data = hulls,
               # Set the transparency of the polygon.
               alpha = 0.15,
               # Set the dimension of the line.
               linewidth = 0.8) +
  # Plot individual observations as semi-transparent circles.
  geom_point(alpha = 0.6,
             size = 2.5) +
  # Plot the centroids into the polygon.
  geom_point(data = centroids, 
             # shape = 15 set the form of the centroids as filled sqaures.
             shape = 15,
             size = 4,
             # Prevents conflicts with the fill aesthetic from the parent ggplot() call.
             inherit.aes = FALSE,
             aes(x = PC1, y = PC2, color = Group)) +
  # Add a label to every point in the plot (with the species name). 
  geom_text_repel(
    aes(label = data$Taxon_name),
    size = 3,
    # Show all the labels without limits.
    max.overlaps = Inf,
    # Set the space around the text.
    box.padding = 0.4,
    # Set the space between the point and its label.
    point.padding = 0.3,
    # Set the colour of the arrow that links point and label.
    segment.color = "grey60",
    # Set the minimum segment lenght.
    min.segment.length = 0.2
  ) +
  # Apply the custom color palette to point borders and polygon outlines.
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  # Apply a black-and-white theme.
  theme_bw() +
  # Label the axes with the variance explained by each component.
  labs(x = paste0("PC1 (", var[1], " %)"), y = paste0("PC2 (", var[2], " %)"))
