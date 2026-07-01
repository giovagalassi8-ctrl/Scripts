#!/usr/bin/env Rscript

# This script calculates the geographic overlap between pairs of species and estimate the amount of geographic space between species within the same clade.
# It construct spatial polygons from occurrence points and computes a Jaccard-like overlap index (intersection area / union area) for every species pair.
# It requires a csv file obtained by the 'geographic_coordinates_gbif.R' script, which contains valid coordinates (WGS84 are preferred).


library(sf)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Import a .csv file obtained by the 'geographic_coordinates_gbif.R' script.
# This file is supposed to contain the following columns: species, decimalLatitude, decimalLongitude.
# (Change with the correct name).
coordinates <- read.csv("Carditidae_coordinates.csv", stringsAsFactors = FALSE)
# Remove invalid coordinates.
coordinates <- coordinates %>%
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude),
         !(decimalLatitude == 0 & decimalLongitude == 0)) %>%
  # Convert to spatial object.
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), 
           # Coordinate Reference System. The 4326 code corresponds to WGS84.
           crs = 4326) %>%
  # Group points by species.
  group_by(species) %>%
  # Count occurrences to filter out species with less than 3 points (cannot form a polygon).
  mutate(n_occurrences = n()) %>%
  filter(n_occurrences >= 3) %>%
  # Create convex hull for each species group and keep the occurrence count.
  summarize(geometry = st_convex_hull(st_union(geometry)),
            n_occurrences = first(n_occurrences)) %>%
  ungroup()

# Alert for species with low occurrences.
# Define the treshold for the minimum number of data points needed (change if necessary).
MIN_OCC <- 5
# Filter the spatial dataframe to isolate species with occurrences below the treshold.
low_occ_sp <- coordinates %>% filter(n_occurrences < MIN_OCC) %>% pull(species)
# Check if the vector contains at least one species name and then print the list of species names that triggered the warning.
if (length(low_occ_sp) > 0) {
  cat("\nWARNING - species with fewer than", MIN_OCC, "occurrences (overlap unreliable):\n")
  print(low_occ_sp)
}

# Extract valid species names.
sp_names <- coordinates$species

# Generate a list of all unique species pairs.
pairs <- combn(sp_names, 2, simplify = FALSE)

# Iterate over pairs to calculate the overlap index using bind_rows for output.
overlap_table <- bind_rows(lapply(pairs, function(p) {
  # Extract geometries for the two species.
  h1 <- coordinates$geometry[coordinates$species == p[1]]
  h2 <- coordinates$geometry[coordinates$species == p[2]]
  
  # Calculate intersection area (suppress topology warnings).
  intersection_area <- suppressWarnings(
    tryCatch(as.numeric(st_area(st_intersection(h1, h2))),
             error = function(e) 0))
  if (length(intersection_area) == 0) intersection_area <- 0
  # Calculate union area.
  union_area <- suppressWarnings(
    tryCatch(as.numeric(st_area(st_union(h1, h2))),
             error = function(e) NA))
  # Return row dataframe.
  data.frame(sp1 = p[1],
             sp2 = p[2],
             overlap_index = intersection_area / union_area)
}))

# Save the final dataframe (change as necessary).
write.csv(overlap_table, "Carditidae_pairwise_overlap.csv", row.names = FALSE)
# Synthesize clade-level overlap.
clade_overlap_mean <- mean(overlap_table$overlap_index, na.rm = TRUE)
cat("\nMean overlap:", round(clade_overlap_mean, 3), "\n")


# Create the final plot to visualize the overlap of the geographic distributions.
get_label_point <- function(geom) {
  # Try to get a point guaranteed to lie on the polygon's surface.
  pt <- tryCatch(st_point_on_surface(geom), error = function(e) NULL)
  # Fallback if the above failed (degenerate/invalid geometry).
  if (is.null(pt) || length(pt) == 0) {
    # Get the bounding box of the polygon.
    bb <- st_bbox(geom)
    # Use the bounding box center as the label anchor point instead.
    pt <- st_sfc(st_point(c(mean(bb[c("xmin", "xmax")]), mean(bb[c("ymin", "ymax")]))),
                 crs = st_crs(geom))
  }
  # Return the anchor point as plain [X, Y] coordinates.
  st_coordinates(pt)
}
# Compute the label anchor point for every species and stack results into one matrix.
label_coords <- do.call(rbind, lapply(coordinates$geometry, get_label_point))
# Store the anchor point coordinates as new columns.
coordinates$label_x <- label_coords[, "X"]
coordinates$label_y <- label_coords[, "Y"]

# Build the plot.
final_plot <- ggplot(coordinates) +
  # Draw each species' range polygon, colored by species and semi-transparent.
  geom_sf(
    aes(fill = species,
        color = species),
    alpha = 0.3,
    linewidth = 0.6) +
  # geom_text_repel automatically nudges overlapping labels apart and draws a thin connecting segment back to the anchor point whenever a label had to be moved away from it.
  geom_text_repel(aes(x = label_x, y = label_y, label = species, color = species),
                  size = 2.2,                # label font size
                  fontface = "italic",       # italics, standard for species names
                  segment.color = "grey40",  # leader line color
                  segment.size = 0.4,        # leader line thickness
                  min.segment.length = 0,    # always draw the connecting segment
                  box.padding = 0.3,         # min spacing kept around each label
                  max.overlaps = Inf,        # never hide a label for being too crowded
                  show.legend = FALSE) +     # no separate legend for this layer
  # Set the plot title and drop the default axis titles.
  labs(title = "Geographic Ranges", x = NULL, y = NULL) +
  theme_minimal() +
  # Remove the color legend since species names are already shown as text.
  theme(legend.position = "none")

# Show the plot.
print(final_plot)
