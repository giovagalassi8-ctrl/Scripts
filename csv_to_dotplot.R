#!/usr/bin/env Rscript

# This script reads a generic tab-delimited text file containing summary statistics.

# Object names are parsed assuming the following pattern:
#     PREFIX_A_PREFIX_B_SUFFIX.out
# Example:
#     concatenated_MS90_g80_allgenes.out
#
#   From this pattern the script extracts:
#     - Group (Y axis): MS90_g80
#     - Color group   : allgenes
#
#   Objects sharing the same group (MS90_g80) are plotted
#   on the same Y-axis line with different colors.

# USAGE:
# [bash] Rscript interactive_grouped_dotplot.R
# [RStudio] source("interactive_grouped_dotplot.R")


library(tidyverse)

# Read the file with the statistics (change with the correct file path).
data <- read_table("<INPUT_FILE>")

# Rename first column as Object for clarity (it's the column containing the name of the objects).
colnames(data_raw)[1] <- "Objects"

# Creates the vector containing the colors to be assigned to the various groups recognized by the script based on the object name
color <- c(
  "allgenes" = "steelblue",
  "rcv" = "tomato",
  "lb" = "forestgreen"
)

# The next stage is data parsing [Expected structure (example): concatenated_MS90_g80_allgenes.out].
data_parsed <- data %>% mutate( 
  
    # Removes the file extension, in this case .out; if it is different, change it accordingly.
    Object_clean = gsub("\\.out$", "", Objects),
    
    # Extract group for Y axis: after first "_" up to third "_"
    Group = sub("^[^_]+_([^_]+_[^_]+)_.*", "\\1", Object_clean),
    
    # Extract suffix for color: after third "_" to end
    ColorGroup = sub(".*_[^_]+_[^_]+_([^_]+)$", "\\1", Object_clean)
  )

# Select only the column needed by the chart, renaming the one of the values (change the name based on the column of interest).
data_plot <- select(data_parsed, Group, ColorGroup, Value = all_of("Missing_percent"))

# Create the dotplot
plot <- ggplot(data_plot, aes(x = Value, y = Group, color = ColorGroup)) +
  # Connects points in the same Group with a black line (shows the distance between values).
  geom_line(aes(group = Group), linetype = "solid", color = "black") +
  # Draw colored dots for ColorGroup.
  geom_point(size = 4, shape = 16) +
  # Apply the colors defined in the color vector.
  scale_color_manual(values = color) +
  # Set the x-axis parameters.
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  # Apply a black-and-white theme to the background.
  theme_bw() +
  # Remove the default grid from the background.
  theme(
    panel.border     = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black")
  ) +
  # Add the labels to the graph
  labs(title = paste("Missing Percent Dotplot"),
       x = "Missing_percent", y = "Group", color = "Type")

