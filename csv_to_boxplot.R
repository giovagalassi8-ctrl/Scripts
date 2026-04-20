#!/usr/bin/env Rscript

# This script reads a CSV file containing different data from different species (each one associated to an higher taxonomic group, e.g. Phylum).
# It produces a boxplot visualizing the distribution of the selected statistic using ggplot2, with custom colors and clean formatting.

# USAGE:
# [Rstudio] source(boxplot_with_table.R)


library(ggplot2)

# Read the CSV file with the interest data (change with the correct .csv file name).
data <- read.csv("INPUT_FILE")

# Rename the first column to "groups" for clarity (it's the column containing the higher taxonomic grouping).
colnames(data)[1] <- "groups"

# Initialize the ggplot object.
ggplot(data = data,aes(
  x = groups, # Taxonomic group on the horizontal axis.
  y = Missing_percent, # Missing data percentage on the vertical axis (change with the name of the column -statistic- of interest).
  fill = groups # Map colours based on groups.
  )) +
  # Add a boxplot layer.
  geom_boxplot(
    outlier.shape = 21, # Set the outlier shape (in this case outliers are shown as filled circles).
    outlier.size  = 2, # Set the oulier point size.
    alpha = 0.8 # Set the transparency (0 -> max, 1-> min).
  ) +
  # Manually assign fill colors to each taxonomic group.
  scale_fill_manual(
    # Each named value maps a group label to a specific color.
    # In the following example, we considered two Gastrotrichs orders and the phylum Platyhelmintes: change accordingly to your data.
    values = c("Chaetonotida" = "blue", "Macrodasyida" = "darkgreen", "Platyhelmintes" = "red")
  ) +
  
  # Set axis and legend labels (change as you desire).
  labs(
    x = "Taxonomic group",
    y = "Missing Percent (%)",
    fill = "groups"
  ) +
  # Apply a black and white theme.
  theme_bw() +
  # Adds more graphic details.
  theme(
    plot.title   = element_text(face = "bold", size = 13, hjust = 0.5),
    axis.text.x  = element_text(size = 11),
    axis.title   = element_text(size = 12),
    legend.position = "none"
  )

