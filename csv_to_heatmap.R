#!/usr/bin/env Rscript

# This script is designed to visualize data from a CSV file as an annotated heatmap. 
# It handles complex layouts by merging multi-row headers into clear column labels and applies a custom Red-White-Blue color gradient (0-50-100) to represent intensity.

# IMPORTANT: PRELIMINARY CHECK OF THE .CSV FILE
# Before running this script, it is recommended to open the .csv file and 
# carefully check the header rows (the first rows defining the columns).
# Sometimes, when converting from Excel to CSV, cells that were "merged" in Excel 
# become one filled cell followed by empty cells in the CSV.
# If there are unwanted white/empty cells in the header rows, the script might 
# generate incomplete labels.

# -> FILL THESE CELLS with the correct text before launching the script so that 
#    every column has a complete and correct label in the plot.

# USAGE:
# [bash] Rscript csv_to_heatmap.R <filename>
# [Rstudio] source(csv_to_heatmap.R)


library(ggplot2)
library(tidyr)
library(dplyr)

# If the file used has more header than 1, we first read the three header rows separately (to assign object's name on the graph).
header1 <- read.csv("01_Heatmap_Matrix/extracted_Table.csv", header = FALSE, nrows = 1)
header2 <- read.csv("01_Heatmap_Matrix/extracted_Table.csv", header = FALSE, skip = 1, nrows = 1)
header3 <- read.csv("01_Heatmap_Matrix/extracted_Table.csv", header = FALSE, skip = 2, nrows = 1)

# Read the data rows (everything after the 3 header rows).
data <- read.csv("01_Heatmap_Matrix/extracted_Table.csv", header = FALSE, skip = 3,
                 colClasses = "character")

# Build x-axis labels by combining the three header rows.
# Remove the first element of each header row, then paste the three rows together with a newline separator.
x_labels <- paste(as.character(row1[-1]),
                  as.character(row2[-1]),
                  as.character(row3[-1]),
                  sep = "\n")

# Extract model names (y-axis) and numeric values.
model_names <- data[, 1]          # first column  → y-axis labels.
values      <- data[, -1]         # all other columns → cell values.

# Convert to a long ("tidy") data frame for ggplot.
# Assign the combined x-axis labels as column names.
colnames(values) <- x_labels

# Add the model names as a column, then reshape from wide to long format.
data_long <- values %>%
  mutate(model = model_names) %>%
  pivot_longer(cols      = -model,
               names_to  = "x_label",
               values_to = "value") %>%
  mutate(
    # Convert values to numbers (empty cells become NA).
    value   = suppressWarnings(as.numeric(value)),
    # Keep the original column order on the x-axis.
    x_label = factor(x_label, levels = x_labels),
    # Keep the original row order on the y-axis (top = first row).
    model   = factor(model,   levels = rev(model_names))
  )

# Draw the heatmap.
heatmap <- ggplot(df, aes(x = x_label, y = model, fill = value)) +
  # Draw tiles; 'color' = border color, 'linewidth' = border thickness.
  # 'width'/'height' < 1 creates the small white gap between tiles.
  geom_tile(color = "black", linewidth = 0.5, width = 0.85, height = 0.85) +
  # Print the numeric value inside each tile (skip NA cells).
  geom_text(aes(label = ifelse(is.na(value), "", value)),
            size = 3, color = "white") +
  # Set the color gradient.
  scale_fill_gradient2(low      = "red",
                       mid      = "white",
                       high     = "blue",
                       midpoint = 50,
                       limits   = c(0, 100),
                       na.value = "grey90",
                       name     = "Value") +
  # Move x-axis labels to the TOP of the plot.
  scale_x_discrete(position = "top") +
  # Axis and plot titles.
  labs(title = "MS90", x = NULL, y = NULL) +
  # Clean white theme.
  theme_void() +
  theme(
    axis.text.x  = element_text(angle = 0, hjust = 0.5,
                                vjust = 0.5, lineheight = 1.1),
    axis.text.y  = element_text(face = "bold"),
    panel.grid   = element_blank(),
    plot.title   = element_text(face = "bold", hjust = 0.5)
  )
