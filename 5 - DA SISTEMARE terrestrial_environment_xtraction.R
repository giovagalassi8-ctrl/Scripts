# =============================================================================
# Script: 06_terrestrial_env_extraction_minimal.R
# =============================================================================

library(terra)
library(geodata)
library(dplyr)

clades <- c("Arion", "Clausiliidae", "Haemadipsidae", "Hormogastridae", "Lumbricidae", "Theba")

# 1. SCARICA E UNISCI I RASTER (Stacking) in un unico oggetto a 3 livelli
env_stack <- c(
  elevation_global(res = 10, path = "wc_data"),
  worldclim_global(var = "bio", res = 10, path = "wc_data")[[c("wc2.1_10m_bio_1", "wc2.1_10m_bio_12")]]
)
names(env_stack) <- c("Elevation", "Temperature", "Precipitation")

# 2. LEGGI TUTTI I CSV INSIEME (Creando una singola grande tabella)
all_coords <- bind_rows(lapply(clades, function(c) {
  file <- paste0(c, "_coordinates.csv")
  if(file.exists(file)) mutate(read.csv(file, stringsAsFactors = FALSE), Clade = c)
})) %>% filter(!is.na(decimalLongitude), !is.na(decimalLatitude))

# 3. ESTRAZIONE SPAZIALE (Una sola operazione vettorializzata per tutti i punti)
points <- vect(all_coords, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
env_data <- extract(env_stack, points)

# 4. CALCOLO MEDIE E SALVATAGGIO
final_summary <- cbind(all_coords, env_data) %>%
  group_by(species, Clade) %>%
  summarize(
    Mean_Elevation = mean(Elevation, na.rm = TRUE),
    Mean_Temperature = mean(Temperature, na.rm = TRUE),
    Mean_Precipitation = mean(Precipitation, na.rm = TRUE),
    N_Occurrences_Env = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(Mean_Elevation)) # Scarta i punti caduti in mare (NA)

write.csv(final_summary, "Phase2_Terrestrial_Environment.csv", row.names = FALSE)
cat("Estrazione completata. File salvato!\n")
