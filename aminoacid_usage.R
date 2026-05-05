#!/usr/bin/env Rscript

# This script generates a heatmap using the amino acid frequencies (previously normalized) of the species under study.
# Amino acids are grouped based on the Dayhoff6 classification, but you can set a random order or an arbitrary one.
# Species (row) are hierarchically clustered using the Ward.D2 method on Euclidean distances. 
# A color sidebar annotates each species by its taxonomic group.

# IMPORTANT:
# Make sure the data is normalized before launching the script; 
# in this case the various amino acid frequencies have been normalized using the following command for each amino acid (columns):
# data <- transform(data, 'column name' = ('column name' / Sum_Aminoacid))
# In particular, we divide each value for the total amino acid usage of the specie (sum of every value on the same row: so, the sum of every row's value equals 1).

# USAGE:
# [Rstudio] source(aminoacid_usage.R)


library(pheatmap)
library(RColorBrewer)
library(viridis)
library(paletteer)

# Read the CSV file into a data frame (change with the correct file name).
data <- read.csv("INPUT_FILE", check.names = FALSE)  # check.names = FALSE prevents R from automatically altering column names

# Rename the first column to "Groups" for clarity.
colnames(data)[1] <- "Groups"

# Create the vector of the 20 amino acid (one-letter code).
# These will be used to select the relevant columns from the data frame.
aa <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

# Define the taxonomic groups of interest.
groups <- c("Platyhelmintes", "Macrodasyida", "Chaetonotida")

# Removes all species in the dataset that are not of interest, 
# filtering them through the groups vector (containing the interest groups instead).
filtered_data <- data[data$Groups %in% groups, ]
# Sort the species in the table by taxonomic group to which they belong.
filtered_data <- filtered_data[order(filtered_data$Groups), ]

# Extract only the amino acid columns and convert to a numeric matrix.
# pheatmap() requires a matrix (not a data frame) as input.
matrix <- as.matrix(filtered_data[, aa])
# gsub() replaces patterns within character sequences.In this case the "_" has been replaced with a space.
rownames(matrix) <- gsub("_", " ", filtered_data$Taxon_name)

# Define amino acid order following the Dayhoff classification.
# Ignore if you don't want to sort amino acids in this way (you have also to set 'cluster_cols=TRUE' into heatmap script).
dayhoff <- c(
  "A","G","P","S","T", #Aliphatic/Polar
  "C", #Cysteine
  "D","E","N","Q", #Acidic
  "F","W","Y", #Aromatic
  "H","K","R", #Basic
  "I","L","M","V" #Aliphatic/Hydrophobic
)

# Reorder the matrix columns according to the Dayhoff scheme. (You can also change 'dayhoff' with the 'aa' object previously created).
# IF THE PREVIOUS STEP WAS NOT DONE, SKIP THIS COMMAND.
matrix <- matrix[, dayhoff]

# Create a dataframe with the group of each species (use later for the sidebar).
rows_annotation <- data.frame(Groups = filtered_data$Groups, 
                              # row.names ensures the annotation aligns correctly with the heatmap rows.
                              row.names = row.names(matrix))

# Define sidebar colors. Each taxonomic group is assigned a distinct color (change with your considered groups).
colors <- list(Groups = c("Chaetonotida"   = "#BED4E9FF",
                          "Macrodasyida"   = "#19647EFF",
                          "Platyhelmintes" = "#FDB927FF"))

# Create a continuous color gradient for the heatmap cells.
# You can also use a dafault palette from the viridis library, such as 'viridis(100)' or 'magma(100)'.
my_palette <- colorRampPalette(c("white", "black"))
heatmap_scale_colors <- my_palette(100)

# Generate the heatmap
pheatmap(
  mat = matrix,
  color = heatmap_scale_colors,
  
  # Select the agglomeration method to be used.
  # Other options are: "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC), "centroid" (= UPGMC).
  clustering_method = "ward.D2",
  
  cluster_rows = TRUE, # Boolean values determining if rows should be clustered or 'hclust' objects.
  clustering_distance_rows = "euclidean", # Distance metric used for row clustering.
  
  cluster_cols = FALSE, # Boolean values determining if columns should be clustered or 'hclust' objects (Set TRUE if dayhoff grouping has not been considered).
  clustering_distance_cols = "euclidean", # Distance metric used for column clustering. In this case is defined but not applied (cols not clustered).
  
  treeheight_row = 50, # The height of a tree for rows, if these are clustered.
  treeheight_col = 50, # The height of a tree for columns, if these are clustered.
  
  annotation_row = rows_annotation, # Data frame mapping each row to its taxonomic group (sidebar).
  annotation_colors = colors, # Colors for the sidebar annotation
  
  fontsize_row = 8, # dimension for row labels (species names).
  fontsize_col = 10, # dimension for column labels (amino acid codes).
  angle_col = 0, # Angle of the column labels (0, 45, 90, 270, 315).
  italic = TRUE, # Set the row labels font.
  
  main = "Heatmap", # Plot title.
  border_color = "black", # Color of the borders around each heatmap cell.
  legend = TRUE, # Display the color scale legend.
  
  # If necessary, removes the column name (in this case overlaps with the first amino acid in the small size image).
  gaps_col = NULL,
  annotation_names_row = FALSE
)
