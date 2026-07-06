# =============================================================================
# Script: 06_terrestrial_env_extraction_minimal.R
# =============================================================================

library(terra)
library(geodata)
library(dplyr)

# I cladi terrestri che ci interessano
clades_target <- c("Arion", "Clausiliidae", "Haemadipsidae", "Hormogastridae", "Lumbricidae", "Theba")

# 1. RICERCA AUTOMATICA NELLE SOTTOCARTELLE (Il trucco del 'recursive = TRUE')
tutti_i_file <- list.files(pattern = "_coordinates\\.csv$", recursive = TRUE, full.names = TRUE)

# Filtriamo solo i file che appartengono ai cladi terrestri
file_terrestri <- tutti_i_file[basename(tutti_i_file) %in% paste0(clades_target, "_coordinates.csv")]

if(length(file_terrestri) == 0) stop("Nessun file _coordinates.csv trovato nelle sottocartelle per i cladi terrestri!")

cat("Trovati", length(file_terrestri), "file nelle sottocartelle. Inizio l'estrazione...\n")

# 2. SCARICA E UNISCI I RASTER (Stacking) in un unico oggetto a 3 livelli
env_stack <- c(
  elevation_global(res = 10, path = "wc_data"),
  worldclim_global(var = "bio", res = 10, path = "wc_data")[[c("wc2.1_10m_bio_1", "wc2.1_10m_bio_12")]]
)
names(env_stack) <- c("Elevation", "Temperature", "Precipitation")

# 3. LEGGI TUTTI I CSV INSIEME 
all_coords <- bind_rows(lapply(file_terrestri, function(f) {
  # Estrae il nome del clade pulendo il nome del file (es: da "cartella/Arion_coordinates.csv" a "Arion")
  clade_name <- sub("_coordinates\\.csv$", "", basename(f))
  read.csv(f, stringsAsFactors = FALSE) %>% mutate(Clade = clade_name)
})) %>% filter(!is.na(decimalLongitude), !is.na(decimalLatitude))

# 4. ESTRAZIONE SPAZIALE (Una sola operazione vettorializzata)
points <- vect(all_coords, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
env_data <- extract(env_stack, points)

# 5. CALCOLO MEDIE E SALVATAGGIO
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
