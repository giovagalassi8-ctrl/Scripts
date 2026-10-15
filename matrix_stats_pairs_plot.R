#!/usr/bin/env Rscript

# This script creates a pairwise comparison plot (scatterplot matrix) of the summary statistics of different supermatrices, using the GGally package.
# Each panel of the lower triangle shows the relationship between two statistics, with a linear regression line for each filtering category (e.g. allgenes, lb, rcv),
# the diagonal shows the distribution of each statistic as a boxplot, and the upper triangle shows the correlation coefficients.
# It requires a CSV file containing one row per supermatrix, with the matrix name, its filtering category and the statistics of interest
# (which can be obtained with 'AMAS.py summary', see the 'AMAS_matrix_stats.sh' script).
# The final plot is saved both in .pdf and .svg format.


library(GGally)
library(ggplot2)
library(dplyr)
library(svglite)

# Read the CSV file containing the statistics of all the supermatrices (change with the correct file path).
data <- read.csv("02_Matrix_Stats/All_matrices_stats.csv", header = TRUE, stringsAsFactors = FALSE)
# Keep only the supermatrices analysed with a specific approach, in this case those whose name ends with '_PM' (partition model; change accordingly).
data_filtered <- data %>%
  filter(grepl("_PM$", matrix))

# Select the columns (statistics) to be compared in the plot (change with the names of the columns of interest).
cols_to_plot <- c("alignment.lenght", "parsimony.informative.sites",
                  "Proportion.variable.sites", "missing.percent")

# Create the vector containing the colors to be assigned to each filtering category (change accordingly).
my_colors <- c("allgenes" = "skyblue4", "lb" = "grey", "rcv" = "indianred3")

# Define the function that draws the panels of the lower triangle:
# a scatterplot of two statistics with a linear regression line (without confidence interval) for each category.
lowerfun <- function(data, mapping){
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_smooth(method = "lm", se = FALSE, size = 0.8) +
    theme_bw()
}
# Define the function that draws the panels of the diagonal: a boxplot of each statistic for each category.
diagfun <- function(data, mapping){
  ggplot(data = data, mapping = mapping) +
    geom_boxplot(alpha = 0.6) +
    theme_bw()
}

# Create the pairwise comparison plot.
p <- ggpairs(data_filtered,
             columns = cols_to_plot,
             # Color points, lines and boxes based on the category column (change accordingly).
             mapping = aes(color = category, fill = category),
             # Apply the functions defined above to the lower triangle and to the diagonal.
             lower = list(continuous = lowerfun),
             diag = list(continuous = diagfun),
             # Use the legend of the first panel as the legend of the whole plot.
             legend = 1) +
  # Apply the colors defined in the my_colors vector.
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  # Set the size of the axis text and of the panel titles, and place the legend at the bottom.
  theme(axis.text = element_text(size = 7),
        strip.text = element_text(size = 8, face = "bold"),
        legend.position = "bottom")

# Save the plot in .pdf and .svg format (change the file names and dimensions as desired).
ggsave("phylogenomic_plot.pdf", plot = p, width = 10, height = 10, device = "pdf")
ggsave("phylogenomic_plot.svg", plot = p, width = 10, height = 10, device = svglite::svglite)
