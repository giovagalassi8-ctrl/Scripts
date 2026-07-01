#!/usr/bin/env Rscript

# This script calculates the geographic overlap between pairs of species and estimate the amount of geographic space between species within the same clade.
# It construct spatial polygons from occurrence points and computes a Jaccard-like overlap index (intersection area / union area) for every species pair.
# It requires a csv file obtained by the 'geographic_coordinates_gbif.R' script, which contains valid coordinates (WGS84 are preferred).


library(sf)
library(dplyr)
library(ggplot2)

# Import a .csv file obtained by the 'geographic_coordinates_gbif.R' script.
# This file is supposed to contain the following columns: species, decimalLatitude, decimalLongitude.
# (Change with the correct name).
coordinates <- read.csv("COORDINATES.CSV", stringsAsFactors = FALSE)
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
  if (length(inter_area) == 0) inter_area <- 0
  # Calculate union area.
  union_area <- suppressWarnings(
    tryCatch(as.numeric(st_area(st_union(h1, h2))),
             error = function(e) NA))
  # Return row dataframe.
  data.frame(sp1 = p[1],
             sp2 = p[2],
             overlap_index = inter_area / union_area)
}))

# Save the final dataframe (change as necessary).
write.csv(overlap_table, "PAIRWISE_OVERLAP.CSV", row.names = FALSE)
# Synthesize clade-level overlap.
clade_overlap_mean <- mean(overlap_table$overlap_index, na.rm = TRUE)
cat("\nMean overlap:", round(clade_overlap_mean, 3), "\n")

# Create the final plot to visualize the overlap of the geographic distributions.
plot(st_geometry(coordinates),
     # Assign an unique fill color to each species (sequence from 1 to the number of total rows).
     col = adjustcolor(1:nrow(coordinates),
                       # Makes the polygons semi-transparent to reveal overlapping areas.
                       alpha.f = 0.3),
     # Assign a unique border to each polygon, matching its fill color.
     border = 1:nrow(coordinates), 
     # Set the title of the graph (change as desired).
     main = "Geographic Ranges")
