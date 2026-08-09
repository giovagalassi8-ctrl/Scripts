# !/usr/bin/env Rscript

# This script extracts terrestrial environmental variables (in particular Elevation, Temperature, and Precipitation) for a predefined list of animal clades.
# It searches for coordinate files recursively, stacks global climate rasters, extracts values for each occurrence, and computes species-level averages.
# It requires the csv files obtained by the 'geographic_coordinates_gbif.R' script. 
# In this case only the terrestrial taxa distribution coordinates has to be considered.


library(terra)
library(geodata)
library(dplyr)

# Define the target terrestrial clades to be analyzed (change with the correct names).
target_clades <- c("Arion", "Clausiliidae", "Haemadipsidae", "Hormogastridae", "Lumbricidae", "Theba")

# Automatically search for coordinate CSV files across all subdirectories.
# The supposed directory structure is: Main directory -> Folders (Phyla) -> Subfolders (Family and Genus) containing a csv coordinates file each.
all_files <- list.files(pattern = "_coordinates\\.csv$", recursive = TRUE, full.names = TRUE)
# Filter the discovered files to include only those matching our target clades.
terrestrial_files <- all_files[basename(all_files) %in% paste0(target_clades, "_coordinates.csv")]

# Halt execution if no matching files are found.
if(length(terrestrial_files) == 0) {
  stop("No _coordinates.csv files found in subdirectories for the target terrestrial clades!")
}

# Download and stack required global raster layers into a single 3-layer SpatRaster utilizing WorldClim data at 10-minute resolution.
env_stack <- c(
  # Fetch the global elevation model (altitude data).
  elevation_global(
    # Set the resolution to 10 minutes of degree (about 340 km² at the equator per pixel).
    # This is a good compromise between file lightness and accuracy for large-scale analysis.
    res = 10,
    # Save the output in a folder named "wc_data" (change as desired).
    path = "wc_data"),
  # Fetch global bioclimatic variables.
  worldclim_global(
    # with "bio" it download a set of 19 standard climate variables.
    var = "bio",
    # Set the resolution to 10 minutes of degree.
    res = 10,
     # Save the output in a folder named "wc_data" (change as desired).
    path = "wc_data"
  # Subset only two variables: Annual Mean Temperature and Annual Precipitation.
  )[[c("wc2.1_10m_bio_1", "wc2.1_10m_bio_12")]]
)
# Rename the raster layers for easier downstream reference.
names(env_stack) <- c("Elevation", "Temperature", "Precipitation")

# Read all CSV files, append the corresponding clade name, and bind them into a single dataframe.
all_coords <- bind_rows(lapply(terrestrial_files, function(f) {
  # Extract the clade name from the filename (e.g., "folder/Arion_coordinates.csv" -> "Arion")
  clade_name <- sub("_coordinates\\.csv$", "", basename(f))
  read.csv(f, stringsAsFactors = FALSE) %>% mutate(Clade = clade_name)
})) %>% 
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude))  # Remove records missing coordinates.

# Convert the compiled coordinate dataframe into a spatial vector object.
points <- vect(all_coords, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")

# Extract environmental raster values for all points in a single vectorized operation.
env_data <- extract(env_stack, points)


# Combine original coordinates with extracted environmental data, calculate means, and export
final_summary <- cbind(all_coords, env_data) %>%  
  # Group the dataset by both 'species' and 'Clade' so that all subsequent calculations are performed for each unique species within its respective clade.
  group_by(species, Clade) %>%
  # Calculate summary statistics for each group.
  summarize(
    # Compute the average elevation, ignoring missing values (NA).
    Mean_Elevation = mean(Elevation, na.rm = TRUE),
    # Compute the average temperature, ignoring missing values (NA).
    Mean_Temperature = mean(Temperature, na.rm = TRUE),
    # Compute the average precipitation, ignoring missing values (NA).
    Mean_Precipitation = mean(Precipitation, na.rm = TRUE),
    # Count the total number of occurrence records for the current species.
    N_Occurrences_Env = n(),
    # Ungroup the resulting dataframe after summarization to prevent unexpected grouping behaviors in any downstream analysis.
    .groups = "drop"
  ) %>% 
  # Retain only the rows where 'Mean_Elevation' is not NA (effectively discarding marine species or points that fell completely outside the raster boundaries).
  filter(!is.na(Mean_Elevation))

# Save the final aggregated data to a CSV file.
output_filename <- "Terrestrial_Environment.csv"
write.csv(final_summary, output_filename, row.names = FALSE)
