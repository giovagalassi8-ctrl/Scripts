library(ggplot2)
library(dplyr)

#Import data
data <- read.csv("all_results_final2.csv", header = TRUE)

#Filters the data while keeping only the columns of interest :
# -Name of the trees,
# -Taxonomy,
# -Identification of the monophyly,
# -Support
df <- data[,c(1,2,5,6)]

# Add a column and create the simplified labels for the final graph.
# In this case it is expected the following structure: "Prefix_NAME.treefile".
df <- df |>
  mutate(tree_label = source_tree |>
           # Remove the Prefix (change accordingly).
           gsub("ML_MS80_", "", x = _) |>
           # Remove the extension (change accordingly).
           gsub("\\.treefile", "", x = _))

clade <- c("Lophotrochozoa","Gnathifera","Rouphozoa","Ecdysozoa","Trochozoa","Chaetognathifera","Tetraneuralia","Kryptrochozoa")

other_clade <- sort(unique(df$Grouping))

order <- c(clade, setdiff(other_clade, clade))
  
# Set the order in which the labels will appear on the axis.
df <- df |>
  mutate(
    tree_label = factor(tree_label, 
                        levels = rev(unique(tree_label)  # Set the elements on the x-axis on alphabetical order.
                                     )),
    Grouping   = factor(Grouping,
                        levels = order)  # Set the elements on the y-axis based on the order previously defined (if you want an alphabetical order, change with "sort(unique(Grouping)").
                                      )

# Create the plot.
p <- ggplot(df, aes(
  x = tree_label,
  y = Grouping)
  )+
  
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
                          Category == "Polyphyletic"),
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
  # Set the axis elements.
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank(),
  )
