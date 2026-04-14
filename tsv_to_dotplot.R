#!/usr/bin/env Rscript

# This script reads a tab-delimited text file containing summary statistics of the concatenated files used in the analysis, which can be obtained with 'AMAS.py summary'.
# Once the statistic of interest has been chosen (corresponding to a column of the file) it generates a dot-plot in which each point corresponds to a different concatenated file.

# Object names are parsed assuming the following pattern:
#     PrefixA_PrefixB_Suffix.out
# Example:
#     concatenated_MS90_g80_allgenes.out
#     PrefixA = concatenated
#     PrefixB = MS90_g80
#     Suffix = allgenes
#
#   From this pattern the script extracts:
#     - Group (Y axis): MS90_g80
#     - Color group   : allgenes
#
#   Objects sharing the same group (MS90_g80) are plotted
#   on the same Y-axis line with different colors.

# USAGE:
# [RStudio] source("interactive_grouped_dotplot.R")


library(ggplot2)
library(tidyverse)

# Read the file with the statistics (change with the correct file path).
data <- read_table("INPUT_FILE")

# Select the column of interest (change with the correct column name).
selected_column <- "COLUMN_NAME"

# Rename first column as Object for clarity (it's the column containing the name of the concatenated files -or any other object considered-).
# Skip this command if it already has a name.
colnames(data)[1] <- "Objects"

# Creates the vector containing the colors to be assigned to the various groups (suffix) recognized by the script based on the object name.
color <- c(
  "allgenes" = "steelblue",
  "rcv" = "tomato",
  "lb" = "forestgreen"
)

# The next stage is data parsing.
data_parsed <- data %>% mutate( 
  
    # Removes the file extension, in this case .out; if it is different, change it accordingly.
    Object_clean = gsub("\\.out$", "", Objects),
    
    # Extract group (PrefixB) for Y axis.
    # In this case the group (MS90_g80) is after first "_" up to third "_" .
    Group = sub("^[^_]+_([^_]+_[^_]+)_.*", "\\1", Object_clean),
    
    # Extract Suffix for color.
    # In this case is after third "_" to end.
    ColorGroup = sub(".*_[^_]+_[^_]+_([^_]+)$", "\\1", Object_clean)
  )

# Select only the column needed for the plot, renaming the one of the values (change the name based on the column of interest).
data_plot <- select(data_parsed, Group, ColorGroup, Value = all_of(selected_column))

# Create the dotplot
plot <- ggplot(data_plot, aes(x = Value, y = Group, color = ColorGroup)) +
  # Connects points in the same Group with a black line (shows the distance between values).
  geom_line(
    aes(group = Group),
    linetype = "solid",
    color = "black") +
  # Draw colored dots for ColorGroup.
  geom_point(
    size = 4,
    shape = 16) +
  # Apply the colors defined in the color vector.
  scale_color_manual(values = color) +
  # Set the x-axis parameters.
  scale_x_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05))) +
  # Apply a black-and-white theme to the background.
  theme_bw() +
  # Remove the default grid from the background.
  theme(
    panel.border     = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black")
  ) +
  # Add the labels to the graph (change the title and axis name).
  labs(
    title = paste("TITLE"),
    x = "X-AXIS_NAME",
    y = "Y-AXIS_NAME",
    color = "Type") +
  # Set the title in the middle.
  theme(plot.title = element_text(hjust = 0.5))
