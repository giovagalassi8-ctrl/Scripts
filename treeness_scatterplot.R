# !/usr/bin/env/ Rscript

# This script allows you to create a scatterplot to graphically represent the results of 'phykit toverr' (a text file containing the values of treeness/RCV, treeness, and RCV for each object).
# In particular, the treeness (measure of the proportion of the total tree length (sum of all branch lengths) that is found on internal branches) is represented on the x-axis, 
# while the Relative Composition Variability (RCV, measure of the average variability in sequence composition among taxa in a sequence alignment) is represented on the y-axis.
# Higher values of treeness/RCV (high treeness, low RCV) are desirable, as they indicate that a gene is likely to be less susceptible to systematic biases.
# Considering this ratio it is possible to see which is the model that presents the least bias: it can be useful for evaluating which model is best among the various ones used to create trees.

# USAGE:
# [Rstudio] source("treeness_scatterplot.R")

library(ggplot2)

# Read the file containing the treeness and RCV values (change with the correct file path).
data <- read.table('FILE_PATH', header = FALSE)
# IF NECESSARY: remove the column containing the treeness/RCV value (usually the second one), which will not be represented on the graph.
data <- data[,-2]

# Change column names for clarity.
colnames(data)[2] <- "treeness"
colnames(data)[3] <- "RCV"

# Create a vector containing the objects belonging to the model you want to highlight on the graph.
# In this case the vector contains the objects belonging to the model  assumed to be the one with with the best treeness/RCV ratio (change if necessary).
kpic <- data[c(49:57),]

# Assigns a color to the selected objects in the previous vector to highlight them in the graph.
color <- ifelse(data$V1 %in% kpic$V1, "red", "black")

# Create the scatterplot.
p <- ggplot(data, aes(
    x=treeness, # Treeness value on the x-axis.
    y=RCV # RCV value on the y-axis.
    ))+
  # Draw the points on the graph.
  geom_point(
    size = 2.5,
    alpha = 0.85,
    col = color
    ) +
  # Add a bisector to the graph.
  geom_abline(
    intercept = 0,
    slope = 1,
    ) +
  # Fix the aspect ratio so that 1 unit on the x-axis equals 1 unit on the y-axis.
  coord_fixed(ratio = 1) +
  # Add a black and white theme. 
  theme_bw() +
  # Set axis limits (in this case they were like this to correctly represent the bisector. Adjust as you want).
  scale_x_continuous(limits = c(0.215, 0.255)) +  
  scale_y_continuous(limits = c(0.215, 0.255)) +
  # Add labels to the graph.
  labs(
    x = "Treeness",
    y = "Relative Composition Variability (RCV)",
    col = NULL
  ) +
  # Remove background grid and add a legend.
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    axis.text = element_text(color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 12)
  )
