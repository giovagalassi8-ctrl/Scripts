#!/bin/env Rscript

# This script reads a CSV file containing the monophyly status and the support value for each taxonomic group of interest within the phylogenetic trees considered,
# and creates an heatmap that allows you to visualize the monophyly of each group. On the x-axis of the graph every phylogenetic trees will be reported, while on
# the y-axis there will be all taxonomic groups.
# I suggest for the imput file to use MonoPhylo.py, that creates a tsv file with all the necessary parameters, but every script that returns the same output is fine:
# it's just important that this file has this expected columns: Tree Name, Taxonomic Group, Monophyly Status and Support Value.


library(ggplot2)
library(dplyr)

# Read the CSV file using MonoPhylo.py with the interest data (change with the correct .csv file name).
data <- read.csv("INPUT_FILE", header = TRUE)

# Filters the data while keeping only the columns of interest :
# -Tree Name,
# -Taxonomic Group,
# -Monophyly Status,
# -Support Value.
# Change with the correct column number.
df <- data[,c(1,2,5,6)]

# Add a column and create the simplified labels for the final graph.
# In this case it is expected the following structure: "Prefix_NAME_pruned.treefile".
df <- df |>
  mutate(tree_label = source_tree |>
           # Remove the Prefix (change accordingly, e.g. "ML_MS90_" for the MS90 supermatrix).
           gsub("ML_MS80_", "", x = _) |>
           # Remove the extension (change accordingly).
           # "[._]pruned" also removes the "_pruned" (or ".pruned") suffix of the pruned trees.
           gsub("[._]pruned\\.treefile", "", x = _))

# Use this line ONLY if you want to remove all the lines where Category has NA across ALL trees in the graph.
df <- df |>
  group_by(Grouping) |>
  filter(!all(is.na(Category))) |>
  ungroup()

# Creates a vector to order the taxonomic groups on the y-axis.
# Change with your groups or skip to use a random oreder.
# The vector is split into three blocks (phyla, groups joined with "+", higher clades),
# so that each block can be shown as a separate panel, divided by a white space.
phyla_vector <- c(
  "Acanthocephala", "Annelida", "Arthropoda", "Brachiopoda", "Bryozoa",
  "Chaetognatha", "Cycliophora", "Entoprocta", "Gastrotricha", "Gnathostomulida",
  "Kinorhyncha", "Micrognathozoa", "Mollusca", "Nematoda", "Nematomorpha",
  "Nemertea", "Onychophora", "Phoronida", "Platyhelminthes", "Priapulida",
  "Rotifera", "Tardigrada")

combined_vector <- c(
  "Annelida+Brachiopoda+Phoronida", "Annelida+Brachiopoda+Phoronida+Mollusca",
  "Annelida+Nemertea", "Annelida+Nemertea+Brachiopoda+Phoronida",
  "Brachiopoda+Phoronida", "Entoprocta+Cycliophora", "Mollusca+Brachiopoda+Phoronida",
  "Mollusca+Nemertea", "Mollusca+Brachiopoda+Phoronida+Nemertea",
  "Platyzoa+Entoprocta+Cycliophora", "Platyzoa+Polyzoa",
  "Trochozoa+Bryozoa", "Trochozoa+Entoprocta+Cycliophora", "Trochozoa+Polyzoa")

clades_vector <- c(
  "Eutrochozoa", "Lophophorata", "Platytrochozoa", "Platyzoa", "Polyzoa",
  "Kryptrochozoa", "Tetraneuralia", "Chaetognathifera", "Trochozoa", "Syndermata",
  "Ecdysozoa", "Rouphozoa", "Gnathifera", "Lophotrochozoa")

taxa_vector <- c(phyla_vector, combined_vector, clades_vector)

# If you want to order the x-axis based on the number of monophyletic groups per matrix,
# create the following object.
mono_order <- df |>
  filter(Category == "Monophyletic") |>
  # Avoid duplicates.
  distinct(tree_label, Grouping) |>
  count(tree_label, name = "n_mono") |>
  arrange(n_mono)

# Set the order in which the labels will appear on the axis.
df <- df |>
  mutate(
    # Set the elements on the x-axis based on the defined classification (in this case: allgenes, rcv and lb).
    block = case_when(
      grepl("allgenes", tree_label) ~ 1,
      grepl("rcv", tree_label) ~ 2,
      grepl("lb", tree_label) ~ 3,
      TRUE ~ 4
    ),
    # Set the elements on the x-axis in alphabetical order accordingly to the previously created groups.
    # If you don't want this grouping, change with 'tree_label = factor(tree_label, levels = rev(unique(tree_label)' and ignore the 'block =' line.
    # If you want to order the matrices labels based on the number of monophyletic groups, change with 'tree_label = factor(tree_label, levels = mono_order$tree_label),' the next command line.
    tree_label = factor(tree_label, levels = unique(tree_label[order(block, tree_label)])),
    
    # Set the elements on the y-axis based on the order previously defined.
    Grouping = factor(Grouping,
                      levels = taxa_vector),
    
    # Add a column with the phylogenetic approach of each tree, read from the end of its label.
    # It defines the vertical panels (one per approach): change the order of 'levels' to change the order of the panels.
    model = factor(case_when(
      grepl("_MM$", tree_label) ~ "MM",
      grepl("_6aa$", tree_label) ~ "6aa",
      grepl("_PM$", tree_label) ~ "PM"),
      levels = c("MM", "6aa", "PM")),
    
    # Add a column with the block of each taxonomic group (phyla, "+" groups, higher clades).
    # It defines the horizontal panels.
    row_block = factor(case_when(
      Grouping %in% phyla_vector ~ "Phyla",
      Grouping %in% combined_vector ~ "Combined groups",
      Grouping %in% clades_vector ~ "Higher clades"),
      levels = c("Phyla", "Combined groups", "Higher clades")),
    
    # Add a column with the class of each cell for the 3-colour plot (2.); the cells with NA are "Not tested".
    # The threshold of 95 is the one usually adopted for UltraFast bootstrap support.
    Status = factor(case_when(
      is.na(Category) ~ "Not tested",
      Category == "Monophyletic" & !is.na(Support) & Support >= 95 ~ "Monophyletic (support ≥ 95)",
      Category == "Monophyletic" ~ "Monophyletic (support < 95)",
      TRUE ~ "Non-monophyletic"),
      levels = c("Monophyletic (support ≥ 95)", 
                 "Monophyletic (support < 95)",
                 "Non-monophyletic", "Not tested")))


# Set the panels and the white spaces between them, shared by both plots.
model_names <- c(MM = "Mixture models", `6aa` = "Dayhoff-6 recoding", PM = "Partition models")

panels <- list(
  # Split the heatmap in panels: blocks of groups in rows, phylogenetic approaches in columns.
  facet_grid(rows = vars(row_block), cols = vars(model),
             # Keep in each panel only its own rows and columns, with all the cells of the same size.
             scales = "free", space = "free",
             labeller = labeller(model = model_names)),
  # Remove the suffix of the approach from the x-axis labels, since it is already written above each panel.
  scale_x_discrete(position = "top", labels = function(x) sub("_(MM|6aa|PM)$", "", x)),
  # Set the width of the white spaces between the panels (x = vertical spaces, y = horizontal spaces).
  theme(panel.spacing.x = unit(4, "mm"),
        panel.spacing.y = unit(3, "mm"),
        strip.text.x = element_text(size = 7, face = "bold"),
        # Hide the names of the blocks of groups on the right side of the heatmap.
        strip.text.y = element_blank(),
        strip.placement = "outside")
)


# THERE ARE 2 TYPES OF POSSIBLE PLOTS:

# 1. Create the plot WITH MONOPHYLY CELLS FILLED WITH SUPPORT VALUES.
p <- ggplot(df, aes(
  x = tree_label,
  # This line sets the y-axis by redefining the order of the Grouping levels in reverse.
  y = factor(Grouping,
             levels = rev(levels(factor(Grouping))))
)) +
  
  # Sets a light gray background on all boxes;
  # it will then be kept only in those whose value is NA.
  geom_tile(
    fill = "grey85",  # Colour of the cell.
    color = "white",  # Colour of the border.
    linewidth = 0.3) +
  
  # Sets a dark gray background on the boxes corresponding to paraphyletic or polyphyletic groups.
  geom_tile(data = filter(df, Category %in% c("Paraphyletic", "Polyphyletic")),
            fill = "grey55",  # Colour of the cell.
            color = "white",  # Colour of the border.
            linewidth = 0.3) +
  # Fill these cells with a "X" (to better indicate that they are not monophyletic).
  geom_text(data = filter(df,
                          # If there are also paraphyletic clades, change with "Category %in% c("Paraphyletic", "Polyphyletic"))," .
                          Category %in% c("Paraphyletic", "Polyphyletic")),
            label = "X",
            size = 1.8,
            color = "white") +
  
  # Sets a dark gray background to cells that correspond to monophyletic groups
  geom_tile(data = filter(df, Category == "Monophyletic"),
            fill = "grey25",  # Colour of the cell.
            color = "white",  # Colour of the border.
            linewidth = 0.3) +
  # Fill these cells with the corresponding Support value.
  geom_text(data = filter(df, Category == "Monophyletic"),
            aes(label = Support),  # Add the support value into the cell.
            size = 1,
            color = "white") +
  
  # Add the labels to the axis
  # (they have been removed so as not to cause confusion. If you want, add it).
  labs(x = NULL, y = NULL) +
  # Add a minimal theme to the graph.
  theme_minimal(base_size = 8) +
  # Set the axis elements.
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank(),
  ) +
  # Add the panels and the white spaces.
  panels


# 2. Create the plot WITH 3 COLOURS CORRESPONDING TO THE MONOPHYLY STATUS AND SUPPORT VALUE.
# Set one colour for each class of the 'Status' column created above; the cells not tested stay white.
# Alternative colours: "#63859C" (support >= 95), "#C9D6DF" (support < 95), "#C8796B" (non-monophyletic).
status_colours <- c(
  "Monophyletic (support ≥ 95)" = "#2C6A94",
  "Monophyletic (support < 95)"     = "#9CC3DE",
  "Non-monophyletic"                = "#C8796B",
  "Not tested"                      = "white")

# Create the plot
d <- ggplot(df, aes(
  x = tree_label,
  # This line sets the y-axis by redefining the order of the Grouping levels in reverse.
  y = factor(Grouping,
             levels = rev(levels(factor(Grouping))))
)) +
  
  # Fill each cell with the colour of its class.
  geom_tile(
    aes(fill = Status),
    color = "white",
    linewidth = 0.3
  ) +
  
  # Apply the colours and keep all the classes in the legend.
  scale_fill_manual(
    values = status_colours,
    drop = FALSE,
    name = NULL,
    # Add a thin grey border around the legend keys, so that the white "Not tested" key is visible.
    guide = guide_legend(override.aes = list(color = "grey70"))
  ) +
  
  # Add the labels to the axis.
  # (they have been removed so as not to cause confusion. If you want, add it).
  labs(x = NULL, y = NULL) +
  # Add a minimal theme to the graph.
  theme_minimal(base_size = 8) +
  # Set the axis elements.
  theme(
    axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank(),
    # Place the legend below the heatmap, on one line.
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.key.size = unit(3, "mm")
  ) +
  # Add the panels and the white spaces.
  panels
