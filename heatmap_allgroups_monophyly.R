# !/bin/env Rscript

# This script reads a CSV file containing the monophyly status and the support value for each taxonomic group of interest within the phylogenetic trees considered,
# and creates an heatmap that allows you to visualize the monophyly of each group. On the x-axis of the graph every phylogenetic trees will be reported, while on
# the y-axis there will be all taxonomic groups.
# I suggest for the imput file to use MonoPhylo.py, that creates a tsv file with all the necessary parameters, but every script that returns the same output is fine:
# it's just important that this file has this expected columns: Tree Name, Taxonomic Group, Monophyly Status and Support Value.

# USAGE:
# [Rstudio] source(heatmap_allgroups_monophyly.R)


library(ggplot2)
library(dplyr)

# Read the CSV file with the interest data (change with the correct .csv file name).
data <- read.csv("INPUT_FILE", header = TRUE)

# Filters the data while keeping only the columns of interest :
# -Tree Name,
# -Taxonomic Group,
# -Monophyly Status,
# -Support Value.
# Change with the correct column number.
df <- data[,c(1,2,5,6)]

# Add a column and create the simplified labels for the final graph.
# In this case it is expected the following structure: "Prefix_NAME.treefile".
df <- df |>
  mutate(tree_label = source_tree |>
           # Remove the Prefix (change accordingly).
           gsub("rooted_ML_MS80_", "", x = _) |>
           # Remove the extension (change accordingly).
           gsub("\\.treefile", "", x = _))

# Use this line ONLY if you want to remove all the lines where Category has NA across ALL trees in the graph.
df <- df |>
  group_by(Grouping) |>
  filter(!all(is.na(Category))) |>
  ungroup()

# Creates a vector to order the taxonomic groups on the y-axis.
# Change with your groups or skip to use a random oreder.
taxa_vector <- c(
  "Acanthocephala", "Annelida", "Arthropoda", "Brachiopoda", "Bryozoa", 
  "Chaetognatha", "Cycliophora", "Entoprocta", "Gastrotricha", "Gnathostomulida", 
  "Kinorhyncha", "Micrognathozoa", "Mollusca", "Nematoda", "Nematomorpha", 
  "Nemertea", "Onychophora", "Phoronida", "Platyhelminthes", "Priapulida", 
  "Rotifera", "Tardigrada",
  
  "Annelida+Brachiopoda+Phoronida", "Annelida+Brachiopoda+Phoronida+Mollusca", 
  "Annelida+Nemertea", "Annelida+Nemertea+Brachiopoda+Phoronida", 
  "Brachiopoda+Phoronida", "Entoprocta+Cycliophora", "Mollusca+Brachiopoda+Phoronida", 
  "Mollusca+Nemertea", "Mollusca+Brachiopoda+Phoronida+Nemertea", "Nemertea+Brachiopoda+Phoronida", 
  "Platyzoa+Entoprocta+Cycliophora", "Platyzoa+Polyzoa", 
  "Trochozoa+Bryozoa", "Trochozoa+Entoprocta+Cycliophora", "Trochozoa+Polyzoa",
  
  "Eutrochozoa", "Lophophorata", "Platytrochozoa", "Platyzoa", "Polyzoa", 
  "Kryptrochozoa", "Tetraneuralia", "Chaetognathifera", "Trochozoa", "Syndermata",
  "Ecdysozoa", "Rouphozoa", "Gnathifera", "Lophotrochozoa")

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
      grepl("rcv",      tree_label) ~ 2,
      grepl("lb",       tree_label) ~ 3,
      TRUE                          ~ 4
    ),
    # Set the elements on the x-axis in alphabetical order accordingly to the previously created groups.
    # If you don't want this grouping, change with 'tree_label = factor(tree_label, levels = rev(unique(tree_label)' and ignore the 'block =' line.
    # If you want to order the matrices labels based on the number of monophyletic groups, change with 'tree_label = factor(tree_label, levels = mono_order$tree_label),' the next command line.
    tree_label = factor(tree_label,
                        levels = unique(tree_label[order(block, tree_label)])),
    
    # Set the elements on the y-axis based on the order previously defined.
    Grouping   = factor(Grouping,
                        levels = taxa_vector))


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
  # Makes cells exactly square.
  coord_fixed() +
  # Add a minimal theme to the graph.
  theme_minimal(base_size = 8) +
  scale_x_discrete(position = "top") +
  # Set the axis elements.
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank(),
  )


# 2. Create the plot WITH A COLOR GRADIENT CORRESPONDING TO THE SUPPORT VALUE.
# Ensure the Support column is treated as numeric for the continuous gradient
df$Support <- as.numeric(df$Support)

# Create the plot
d <- ggplot(df, aes(
  x = tree_label,
  # This line sets the y-axis by redefining the order of the Grouping levels in reverse.
   y = factor(Grouping,
             levels = rev(levels(factor(Grouping))))
)) +
  
  # Sets a light gray background on all boxes; 
  # it will then be kept only in those whose value is NA.
  geom_tile(
    fill = "grey85",
    color = "white",
    linewidth = 0.3
  ) +
  
  # Sets a dark gray background on the boxes corresponding to paraphyletic or polyphyletic groups.
  geom_tile(
    data = filter(df, Category %in% c("Paraphyletic", "Polyphyletic")),
    fill = "rosybrown",
    color = "white",
    linewidth = 0.3
  ) +
  
  # Map the fill color to the Support value for Monophyletic groups.
  geom_tile(
    data = filter(df, Category == "Monophyletic"),
    aes(fill = Support),  # Replaces the fixed fill color
    color = "white",
    linewidth = 0.3
  ) +
  
  # Apply a continuous blue color scale (from light to dark).
  scale_fill_gradient(
    low = "lightblue",   # Color for the lowest support values.
    high = "darkblue",   # Color for the highest support values.
    name = "Support",
    na.value = "transparent"  # Ensures it doesn't overwrite the NA background.
  ) +
  
  # Add the labels to the axis. 
  # (they have been removed so as not to cause confusion. If you want, add it).
  labs(x = NULL, y = NULL) +
  # Makes cells exactly square.
  coord_fixed() +
  # Add a minimal theme to the graph.
  theme_minimal(base_size = 8) +
  scale_x_discrete(position = "top") +
  # Set the axis elements.
  theme(
    axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7)
  )
