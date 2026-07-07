# =============================================================================
# Script: 07_marine_env_extraction_minimal.R
# Purpose: Automatically downloads global marine layers (Bathymetry and Bottom 
#          Temperature) from Bio-ORACLE, extracts values for marine species 
#          occurrences, and calculates species-specific environmental means.
# =============================================================================

library(sdmpredictors)
library(terra)
library(dplyr)

# 1. TARGET MARINE CLADES 
clades_target <- c("Carditidae", "Conidae", "Siboglinidae", "Strombidae", 
                   "Terebratulina", "Terebridae", "Veneridae", "Monocelididae", "Lingula")

# 2. FIND COORDINATE FILES IN SUBDIRECTORIES
all_files <- list.files(pattern = "_coordinates\\.csv$", recursive = TRUE, full.names = TRUE)
marine_files <- all_files[basename(all_files) %in% paste0(clades_target, "_coordinates.csv")]

if(length(marine_files) == 0) stop("No _coordinates.csv files found for marine clades!")
cat("Found", length(marine_files), "marine clades. Starting extraction...\n")

# 3. DOWNLOAD BIO-ORACLE MARINE LAYERS
cat("Downloading/Loading Bio-ORACLE marine layers (Bathymetry & Bottom Temp)...\n")
# BO_bathymean = Bathymetry (Depth in meters)
# BO2_tempmean_bdmean = Sea Bottom Temperature (Celsius)
marine_layers <- load_layers(c("BO_bathymean", "BO2_tempmean_bdmean"))

# Convert to 'terra' SpatRaster format for fast extraction
env_stack <- rast(marine_layers)
names(env_stack) <- c("Depth", "Bottom_Temperature")

# 4. LOAD AND BIND ALL CSVs
cat("Merging coordinates...\n")
all_coords <- bind_rows(lapply(marine_files, function(f) {
  clade_name <- sub("_coordinates\\.csv$", "", basename(f))
  read.csv(f, stringsAsFactors = FALSE) %>% mutate(Clade = clade_name)
})) %>% filter(!is.na(decimalLongitude), !is.na(decimalLatitude))

# 5. SPATIAL EXTRACTION 
cat("Extracting environmental values (this may take a moment)...\n")
points <- vect(all_coords, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
env_data <- extract(env_stack, points)

# 6. CALCULATE SPECIES MEANS & CLEAN UP
final_summary <- cbind(all_coords, env_data) %>%
  group_by(species, Clade) %>%
  summarize(
    Mean_Depth = mean(Depth, na.rm = TRUE),
    Mean_Bottom_Temperature = mean(Bottom_Temperature, na.rm = TRUE),
    N_Occurrences_Env = n(),
    .groups = "drop"
  ) %>%
  # IMPORTANT: Filter out NA depths. Coastal GBIF points (like beach-combed shells) 
  # often fall on "land" pixels where the ocean dataset has no data.
  filter(!is.na(Mean_Depth))

# 7. EXPORT RESULTS
write.csv(final_summary, "Phase2_Marine_Environment.csv", row.names = FALSE)
cat("=== EXTRACTION COMPLETE ===\n")
cat("Extracted data for", nrow(final_summary), "unique marine species.\n")
cat("Saved as 'Phase2_Marine_Environment.csv'\n")

