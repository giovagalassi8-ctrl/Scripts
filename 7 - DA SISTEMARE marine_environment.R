#!/usr/bin/env Rscript

# This script extracts marine environmental variables ( in particular Bathymetry/Depth and Sea Bottom Temperature) for a predefined list of marine clades.
# It automatically downloads global raster layers from Bio-ORACLE, extracts values for each species occurrence, and computes species-specific averages.
# It requires the csv files obtained by the 'geographic_coordinates_gbif.R' script. 
# In this case only the marine taxa distribution coordinates has to be considered.


library(sdmpredictors)
library(terra)
library(dplyr)

# Define the target marine clades to be analyzed (change with the correct names).
clades_target <- c("Carditidae", "Conidae", "Siboglinidae", "Strombidae", "Terebratulina", "Terebridae", "Veneridae", "Monocelididae", "Lingula")

# Automatically search for coordinate CSV files across all subdirectories.
# The supposed directory structure is: Main directory -> Folders (Phyla) -> Subfolders (Family and Genus) containing a csv coordinates file each.
all_files <- list.files(pattern = "_coordinates\\.csv$", recursive = TRUE, full.names = TRUE)
# Filter the discovered files to include only those matching our target marine clades.
marine_files <- all_files[basename(all_files) %in% paste0(clades_target, "_coordinates.csv")]

# Halt execution if no matching files are found.
if(length(marine_files) == 0) {
  stop("No _coordinates.csv files found in subdirectories for the target marine clades!")
}

# Download or load specific marine environmental layers from the Bio-ORACLE database.
# BO_bathymean = Bathymetry (Depth in meters).
# BO2_tempmean_bdmean = Sea Bottom Temperature (Celsius).
marine_layers <- load_layers(c("BO_bathymean", "BO2_tempmean_bdmean"))

# Convert the downloaded layers into a 'terra' SpatRaster object for faster extraction.
env_stack <- rast(marine_layers)
# Rename the raster layers for easier downstream reference.
names(env_stack) <- c("Depth", "Bottom_Temperature")

# Read all CSV files, append the corresponding clade name, and bind them into a single dataframe.
all_coords <- bind_rows(lapply(marine_files, function(f) {
  # Extract the clade name from the filename (e.g., "folder/Conidae_coordinates.csv" -> "Conidae").
  clade_name <- sub("_coordinates\\.csv$", "", basename(f))
  read.csv(f, stringsAsFactors = FALSE) %>% mutate(Clade = clade_name)
})) %>% 
  # Remove records missing spatial coordinates.
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude))

# Convert the compiled coordinate dataframe into a spatial vector object.
points <- vect(all_coords, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
# Extract environmental raster values for all points in a single vectorized operation.
env_data <- extract(env_stack, points)

# Combine original coordinates with extracted environmental data and calculate means.
final_summary <- cbind(all_coords, env_data) %>%
  # Group the dataset by both 'species' and 'Clade'.
  group_by(species, Clade) %>%
  # Calculate summary statistics for each species.
  summarize(
    # Calculate the average bathymetric depth for the current species, ignoring missing values (NA).
    Mean_Depth = mean(Depth, na.rm = TRUE),
    # Calculate the average sea bottom temperature across all valid occurrence points, ignoring missing values (NA).
    Mean_Bottom_Temperature = mean(Bottom_Temperature, na.rm = TRUE),
    # Count the total number of valid occurrence records that contributed to these averages.
    # The n() function counts the number of rows in the current group.
    N_Occurrences_Env = n(),
    # Explicitly drop all grouping structures from the resulting dataframe.
    # Dropping all groups prevents silent logic errors in downstream analyses.
    .groups = "drop"
  ) %>%
  
  # IMPORTANT: Filter out NA depths.
  # Coastal GBIF points (like beach-combed shells) often fall on "land" pixels where the ocean dataset has no data.
  filter(!is.na(Mean_Depth))

# Save the final aggregated data to a CSV file.
output_filename <- "Marine_Environment.csv"
write.csv(final_summary, output_filename, row.names = FALSE)
