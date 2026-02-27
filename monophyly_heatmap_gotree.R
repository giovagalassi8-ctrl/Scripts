#!/usr/bin/env Rscript

#This script takes the result of gotree analysis on monophyly (in csv format) to create a costum heatmap where the various squares are separated from each other.

# Load required packages
library(ggplot2)
library(tidyr)

# Import the CSV file
monophyly <- read.csv("03_Gotree/00_stats/monophyly_results_gotree.csv", stringsAsFactors = FALSE)

# The ggplot2 library strictly requires that data be organized in a "long" format, where each column represents a single logical variable.
# In a classic CSV table (wide format), however, information is mixed between row names, column names, and values in cells.
# Using ggplot2 is mandatory to switch from "matrix" to "3-column" shape.
monophyly_long <- pivot_longer(monophyly, cols = -Tree, names_to = "Group", values_to = "Status")

# Generate the heatmap
heatmap_plot <- ggplot(monophyly_long, aes(x = Tree, y = Group, 
                                           #Set the fill color of the tiles based on the 'Status' column, treating it as character.
                                           fill = as.character(Status))) +  
  # Add a layer to draw rectangles (tiles) for each data point.
  geom_tile(color = "black", width = 0.85, height = 0.85) +
  # Add a title to the plot
  labs(title = "MS90")
  # Assign fill colors directly to the words from Status column.
  scale_fill_manual(values = c("false" = "white", "true" = "darkgreen", 
                               "FALSE" = "white", "TRUE" = "darkgreen")) +
  # Costumize x-axis(move the x-axis and its labels to the top of the plot instead of the default bottom).
  scale_x_discrete(position = "top") +
  # Apply a completely empty theme.
  theme_void() + 
  theme(
  # Customize the text of the top x-axis.
  axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5),  # To remove all the names on the x-axis change with 'element_blank()'
  # Customize the text of the y-axis
  axis.text.y = element_text(hjust = 1, margin = margin(r = 5))
  # Center the title relative to the entire plot area
  plot.title.position = "plot",
  # Edit the title  
  plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 1))
  )

# ggsave settings can be corrected to change the plot size as desired.
ggsave(file = "monophyly_heatmap_gotree.pdf", heatmap_plot, width = 7, height = 3, dpi = 300)
