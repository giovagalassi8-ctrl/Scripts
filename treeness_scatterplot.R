# !/usr/bin/env/ Rscript

# This script allows you to create a scatterplot to graphically represent the results of 'phykit toverr' (a text file containing the values of treeness/RCV, treeness, and RCV for each object).
# In particular, the treeness (measure of the proportion of the total tree length (sum of all branch lengths) that is found on internal branches) is represented on the x-axis, 
# while the Relative Composition Variability (RCV, measure of the average variability in sequence composition among taxa in a sequence alignment) is represented on the y-axis.
# Higher values of treeness/RCV (high treeness, low RCV) are desirable, as they indicate that a gene is likely to be less susceptible to systematic biases.
# Considering this ratio it is possible to see which is the model that presents the least bias: it can be useful for evaluating which model is best among the various ones used to create trees.
# It is possible to set some parameters in the script to visualize other statistics and improve resolution.
# In the final part of the script it is possible to add two density plot (or any kind of plot you desire).

# USAGE:
# [Rstudio] source("treeness_scatterplot.R")


library(ggplot2)
library(patchwork)

# Read the file containing the treeness and RCV values (change with the correct file path).
data <- read.table('INPUT_TSV_FILE', header = TRUE)

# Create the scatterplot.
p <- ggplot(data2, aes(
  x=Treeness,  # Treeness value on the x-axis.
  y=RCV,  # RCV value on the y-axis.
  fill = meansupport,  # Color each point based on the selected statistic (change as desired with the column name containing the selected stats).
  # If you don't want the size of the points depending on a selected value, ignore the following 'size' setting,
  # and set following 'geom_point(size)' for a default dimension of all points.
  size = MissingPercent  # Size of each point depends on the selected statistic (change as desired).
))+
  # Draw the points on the graph.
  geom_point(
    #size = 4.5, # Set the default dimension of each point. Ignore if you want the size to scale based on a certain value.
    alpha = 1,
    shape = 21,  # The shape 21 is the one where points have a border.
    colour = "black",  # The color is referred to the border of the points.
    stroke = 1  # Set border width.
  ) +
  # Creates a color gradient to fill the points based on the selected value.
  scale_fill_steps(
    low = "white", 
    high = "darkred",
    breaks = c(94:100),  # Set the number of colors that constitute the gradient (change accordingly).
    labels = c("94", "", "", "", "", "", "100")  # Display only these values on the legend (change accordingly).
  ) +
  # Consider if the point's dimension depends on a selected parameter.
  # Set the point dimension parameters.
  scale_size_continuous(
    # Smaller values make the size larger and vice versa.
    trans = "reverse",  # Ignore if you want that higher value correspond to a bigger point size, and vice versa.
    name = "Missing Percent"  # Set the name of the legend.
  ) +
  # Add four quadrants to the graph.
  geom_vline(xintercept = 0.235) +
  geom_hline(yintercept = 0.235) +
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
    fill = "Mean Support"
  ) +
  # Remove background grid and add a legend.
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "lightgrey"),
    axis.text = element_text(color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 10)
  )


# The following part can be done in case you want to add other graphs, containing other statistics, to the plot.
# Change plot type and statistics considered as you want.

# Create the density plot for the mean support for the final plot.
a <- ggplot(data2, aes(meansupport)) + 
  geom_density(linewidth = 1) +
  theme_minimal() +
  scale_x_continuous(limits = c(90,100)) +
  labs(
    x = "Mean support",
    y = "density"
  )

# Create the density plot for the mean branch lenght for the final plot.
b <- ggplot(data2, aes(meanbrlen)) + 
  geom_density(linewidth = 1) +
  theme_minimal() +
  scale_x_continuous(limits = c(0,0.5)) +
  labs(
    x = "Mean branch length",
    y = "density"
  )

# Merge the two density plot so that they are one above the other.
c <- a/b

# Combines the scatterplot with the two density plot.
# In the final graph the scatterplot will be on the left, while the two density plot are on the right.
final_plot <- (p-c) + 
  # The legends are placed together at the top of the graph, so that they are clearly visible.
  plot_layout(guides = "collect") & theme(legend.position = "top")
