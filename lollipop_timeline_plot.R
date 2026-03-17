#!/usr/bin/env Rscript

# This script creates a timeline bar on which are collocated different points associated to an event (e.g. different pubblications about the same topic).
# The hight of the different segments is based on a selected variable (e.g. the number of phyla considered in that study),
# while the color of each point is assigned using a continuos gradient that consider a different variable (e.g. the number of species used).

# USAGE:
# [Rstudio] source(lollipop_timeline_plot.R)


library(ggplot2)
library(tidyverse)
library(viridis)

# Read the CSV file with the interest data (change with the correct .csv file name).
data <- read.csv("99_Paper_Plots/00_literature/00_data/Literature.csv", header=FALSE, stringsAsFactors = FALSE)

# IF NECESSARY: transpose the dataset for a better manipulation.
transposed_data <- as.data.frame(t(data), stringsAsFactors = FALSE)

# Extract specific rows into individual vectors. Change if the data are on different columns, instead of rows.
# In this case every row has an intestation, which is removed with [ ,-1].
objects <- as.character(unlist(data[1,-1]))
time <- as.numeric(unlist(data[2,-1]))
phyla <- as.numeric(unlist(data[4,-1]))
# The following row is extracted for color mapping.
number_of_species <- as.numeric(unlist(data[3,-1])) 

# Create the main data frame for plotting.
data_plot <- data.frame(
  Objects = objects,
  Time = time,
  Phyla = phyla,
  Color = number_of_species
)

# Create the lollipop plot.
# Map the 'Time' variable to the x-axis and 'Phyla' to the y-axis globally for all layers.
lollipop_plot <- ggplot(data_plot, aes(x = Time, y = Phyla, group = Objects)) +
  geom_segment(
    aes(x = Time, y = 0, yend = Phyla), 
    color = "gray50", linewidth = 1,
    # 'position_dodge' preserves the vertical position of a 'geom' while adjusting the horizontal position.
    # (In this case is useful to divide overlapping lollipops).
    position = position_dodge(width = 0.3)) +
  # Add the points at the (x, y) coordinates.
  geom_point(
    aes(color = Color), size = 4,
    position = position_dodge(width = 0.3)) + 
  # Add text labels.
  geom_text(
    aes(label = Objects), vjust = -1, size = 3, color = "black",
    position = position_dodge(width = 1.3)) +
  # Apply a continuos color scale.
  scale_color_viridis_c(option = "viridis", direction = -1, limits = c(0, 100)) +
  # Plot, axis and legend titles.
  labs(
    title = "Lollipop Plot",
    x = "Timeline",
    y = "Number of Phyla",
    color = "Species"
  ) +
  theme_bw() +
  # Remove the default grid from the background.
  theme(
    panel.border     = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black")
  ) +
  # Set the x-axis limits.
  scale_y_continuous(limits = c(0,30), breaks = seq(0,30, by = 5))
