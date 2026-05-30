# !/usr/bin/env Rscript

# Generates a horizontal stacked bar chart visualizing BUSCO completeness results across multiple samples/assemblies.
# Each bar represents one sample and is divided into four BUSCO categories: Single-copy (S), Duplicated (D), Fragmented (F),
# and Missing (M), each displayed as a percentage of the total expected BUSCOs.
# It requires a TSV file produced by the BUSCO_summary.R script, containing one row per sample. 


library(ggplot2)
library(reshape2)

# Read the TSV file obtained using the BUSCO_summary.R script (change with the correct name).
busco_clade <- read.table("01_all_short_summaries_table_with_clade.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Transform the object into a long format (molten).
# Each original row (one sample) is expanded into 4 rows, one per BUSCO category.
busco_long <- melt(busco_clade,
                   # Defines which columns should remain fixed (identifier variables) while the rest are stacked (change if necessary).
                   id.vars = c("sample", "clade"),
                   # Columns to collapse into a single column (S,D,F and M identiy the BUSCO results categories).
                   measure.vars = c("S","D","F","M"),
                   # Name given to the new column that holds the category labels (S/D/F/M).
                   variable.name = "BUSCO_category",
                   # Name given to the new column that holds the corresponding percentages.
                   value.name = "percentage")

# Change the category of the just created BUSCO_category column into a factor.
busco_long$BUSCO_category <- factor(
    busco_long$BUSCO_category,
    levels = c("M", "F", "D", "S")
)

# Creates specific labels combining the clade column with the sample column. 
busco_long$label <- paste(busco_long$clade, busco_long$sample, sep = " | ")

# Change the category of the labels vector into a factor, and then sort each label by clade.
busco_long$label <- factor(
    busco_long$label,
    levels = unique(busco_long$label[order(busco_long$clade)])
)

# Assign a specific colour to every busco category.
busco_colors <- c(
    "S" = "#1f77b4",  # Single-Copy in blue.
    "D" = "#aec7e8",  # Duplicated in light-blue.
    "F" = "#ffdb58",  # Fragmented in yellow.
    "M" = "#d62728"   # Missing data in red.
)

# Creates the plot.
p <- ggplot(busco_long, aes(x = label,  
                            y = percentage,
                            # Set the fill colours based on the BUSCO_category column (change accordingly).
                            fill = BUSCO_category)) +
  # Creates the barplot.
  geom_bar(
    # The statistical transformation to use on the data for this layer.
    stat = "identity", 
    # Set the width between two bars (in this case we set a low value).
    width = 0.95) +
  # Manually colours the plot using the previously created vector.
  scale_fill_manual(
    values = busco_colors,
    breaks = c("S", "D", "F", "M"),
    labels = c("Single-copy", "Duplicated", "Fragmented", "Missing")
  ) +
  # Set the axis labels.
  labs(x = "", y = "BUSCO (%)", fill = "BUSCO category") +
  # Apply a minimal theme.
  theme_minimal(base_size = 11) +
  theme(
    # Customize the appearance of the y-axis labels.
    axis.text.y = element_text(
      # Font size in points. Kept small to avoid overlapping when many species are listed on the y-axis.
      size = 6,
      # Line spacing multiplier for multi-line labels.
      lineheight = 0.7,
      # In this case, the horizontal justification set right-aligned labels (hjust = 1).
      hjust = 1,
      # Adds 2 pt of space on the RIGHT side of each label, creating a small gap between the text and the axis line.
      margin = margin(r = 2)
    ),
    # Customize the appearance of the x-axis labels.
    axis.text.x = element_text(size = 10),
    # Define the legend characteristics.
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    # Set the margins of the plot.
    plot.margin = margin(10, 40, 10, 10)
  ) +
  # Flip coordinates so that horizontal becomes vertical.
  coord_flip() +
  guides(fill = guide_legend(reverse = FALSE))
