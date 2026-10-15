#!/usr/bin/env Rscript

# This script compares the summary statistics of different supermatrices grouped by filtering category (allgenes, rcv, lb).
# It produces two figures:
#   1. a grid of boxplots (one per statistic), in which every supermatrix is shown as a jittered point, with a single shared legend;
#   2. a scatterplot of Treeness vs RCV (Relative Composition Variability).
# In both figures the color of the points depends on the filtering category, the shape on the supermatrix (MS80 or MS90),
# and the transparency on the phylogenetic approach (Partition model, Mixture model or Dayhoff-6 recoding).
# It requires a CSV file containing one row per supermatrix, with the matrix name, its filtering category and the statistics of interest
# (which can be obtained with 'AMAS.py summary', 'phykit toverr' and 'gotree stats').
# Matrix names are parsed assuming the following pattern:
#     ML_MS{80|90}_{trimming}_{category}_{6aa|MM|PM}
#   where the suffix identifies the approach:
#     6aa = Dayhoff-6 recoding (see 'run_Dayhoff6.sh')
#     MM  = mixture model (see 'run_Mixture_Models.sh')
#     PM  = partition model (see 'run_Partitions_Model.sh')
# The figures are saved in .pdf, .svg and .png format.


library(ggplot2)
library(patchwork)

# Read the CSV file containing the statistics of all the supermatrices (change with the correct file path).
data <- read.csv("02_Matrix_Stats/All_matrices_stats.csv", stringsAsFactors = FALSE)
# Identify the supermatrix (MS80 or MS90) from the matrix name.
data$supermatrix <- ifelse(grepl("MS80", data$matrix), "MS80", "MS90")
# Identify the approach from the suffix of the matrix name.
data$approach <- NA
data$approach[grepl("_6aa$", data$matrix)] <- "Dayhoff-6"
data$approach[grepl("_MM$",  data$matrix)] <- "Mixture"
data$approach[grepl("_PM$",  data$matrix)] <- "Partition"
# Convert the grouping columns into factors, setting the order in which their levels are shown in the plots.
data$category    <- factor(data$category,    levels = c("allgenes", "rcv", "lb"))
data$supermatrix <- factor(data$supermatrix, levels = c("MS80", "MS90"))
data$approach    <- factor(data$approach,    levels = c("Partition", "Mixture", "Dayhoff-6"))

# Create the vector containing the colors to be assigned to each filtering category (change accordingly).
cat_colors <- c("allgenes" = "blue",
                "rcv"      = "red",
                "lb"       = "darkgrey")
# Assign a point shape to each supermatrix (16 = filled circle, 17 = filled triangle).
shape_map <- c("MS80" = 16, "MS90" = 17)
# Assign a transparency level to each approach (0 -> max, 1 -> min).
alpha_map <- c("Partition" = 1, "Mixture" = 1, "Dayhoff-6" = 0.5)
# Define a common theme for all the plots: black and white theme without the minor grid, legend at the bottom and light grey panel titles.
theme_clean <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey95", color = NA))
# Select the columns (statistics) to be plotted (change with the names of the columns of interest).
params <- c("alignment.lenght",
            "parsimony.informative.sites",
            "Proportion.variable.sites",
            "missing.percent",
            "Treeness",
            "RCV",
            "meansupport",
            "sumbrlen")
# Set the y-axis label of each statistic (change accordingly).
param_labels <- c("alignment.lenght" = "Alignment length (bp)",
                  "parsimony.informative.sites" = "Parsimony-informative sites",
                  "Proportion.variable.sites" = "Proportion variable sites",
                  "missing.percent" = "Missing data (%)",
                  "Treeness" = "Treeness",
                  "RCV" = "RCV",
                  "meansupport" = "Mean bootstrap support",
                  "sumbrlen" = "Sum of branch lengths")

# Create one boxplot for each selected statistic.
plot_list <- lapply(params, function(p) {
  ggplot(data, aes(x = category, y = .data[[p]])) +
    # Draw the boxplots, hiding the outliers since all the points are shown by the following layer.
    geom_boxplot(aes(fill = category, color = category),
                 alpha = 0.15, outlier.shape = NA, width = 0.6) +
    # Draw every supermatrix as a jittered point (the seed makes the jitter reproducible).
    geom_jitter(aes(color = category, shape = supermatrix, alpha = approach),
                position = position_jitter(width = 0.18, seed = 42),
                size = 2) +
    # Apply the colors, shapes and transparency levels defined above.
    scale_fill_manual(values = cat_colors,  guide = "none") +
    scale_color_manual(values = cat_colors, guide = "none") +
    scale_shape_manual(values = shape_map) +
    scale_alpha_manual(values = alpha_map) +
    # Set the y-axis label and remove the x-axis one.
    labs(x = NULL, y = param_labels[[p]]) +
    theme_clean +
    # Remove the legend from the single plots (a shared one is added below).
    theme(legend.position = "none")
})
# Create a plot used only to extract the legend shared by all the boxplots.
legend_plot <- ggplot(data, aes(x = category, y = Treeness)) +
  geom_point(aes(color = category, shape = supermatrix, alpha = approach), size = 3) +
  scale_color_manual(values = cat_colors, name = "Filtering scheme") +
  scale_shape_manual(values = shape_map,  name = "Supermatrix") +
  scale_alpha_manual(values = alpha_map,  name = "Approach") +
  theme_clean +
  theme(legend.position = "bottom")

# Extract the shared legend (the cowplot package must be installed).
shared_legend <- cowplot::get_legend(legend_plot)
# Arrange the boxplots in a grid with four columns and add the shared legend at the bottom.
boxplot_grid <- wrap_plots(plot_list, ncol = 4) /
  shared_legend +
  plot_layout(heights = c(20, 1))
# Create the scatterplot of Treeness vs RCV.
scatter_plot <- ggplot(data,
                       aes(x = RCV, y = Treeness,
                           color = category, shape = supermatrix, alpha = approach)) +
  geom_point(size = 3) +
  scale_color_manual(values = cat_colors, name = "Filtering scheme") +
  scale_shape_manual(values = shape_map,  name = "Supermatrix") +
  scale_alpha_manual(values = alpha_map,  name = "Approach") +
  labs(x = "RCV (relative compositional variability)",
       y = "Treeness") +
  theme_clean +
  theme(legend.position = "bottom")

# Set the output directory (change accordingly) and create it if it does not exist.
outdir <- "./plots"
dir.create(outdir, showWarnings = FALSE)

# Save both figures in .pdf, .svg and .png format (change dimensions and resolution as desired).
for (ext in c("pdf", "svg", "png")) {
  ggsave(file.path(outdir, paste0("matrix_boxplots.", ext)),
         boxplot_grid, width = 13, height = 7, dpi = 300)
  ggsave(file.path(outdir, paste0("treeness_vs_rcv.", ext)),
         scatter_plot, width = 13, height = 7, dpi = 300)
}
