# Installa i pacchetti se non li hai già: install.packages(c("ape", "sf", "dplyr"))
library(ape)      # Per manipolare l'albero filogenetico
library(sf)       # Per la gestione dei dati spaziali e la creazione dei poligoni
library(dplyr)    # Per manipolare il database in modo semplice

# 1. Definiamo i file (assicurati che siano nella tua working directory di R)
tree_file <- "Carditida.nwk"
gbif_file <- "Carditida_coordinates.csv"

# 2. Carichiamo l'albero filogenetico
phylo_tree <- read.tree(tree_file)

# 3. Carichiamo il dataset spaziale pulito (ora possiamo usare il read.csv normale!)
occurrences <- read.csv(gbif_file, stringsAsFactors = FALSE)

# 4. Filtriamo i dati spaziali e formattiamo i nomi
cleaned_occurrences <- occurrences %>%
  # Teniamo solo le righe che hanno coordinate valide
  filter(!is.na(decimalLatitude) & !is.na(decimalLongitude)) %>%
  # Teniamo solo quelle identificate a livello di specie
  filter(taxonRank == "SPECIES") %>%
  # Sostituiamo gli spazi con l'underscore per farli combaciare con l'albero!
  mutate(formatted_species = gsub(" ", "_", species))

# 5. Troviamo l'intersezione (quali specie abbiamo in entrambi i file?)
gbif_species_list <- unique(cleaned_occurrences$formatted_species)
tree_species_list <- phylo_tree$tip.label
shared_species <- intersect(tree_species_list, gbif_species_list)

cat("--- RISULTATI DELL'INTERSEZIONE ---\n")
cat("Specie filogenetiche originali:", length(tree_species_list), "\n")
cat("Specie con coordinate valide:", length(gbif_species_list), "\n")
cat("Specie condivise per l'analisi (Overlap):", length(shared_species), "\n\n")

# 6. "Pruning": Tagliamo l'albero per mantenere solo le specie condivise
pruned_tree <- keep.tip(phylo_tree, shared_species)
cat("Albero filogenetico tagliato con successo.\n")

# 7. Modulo Spaziale: Creiamo i Poligoni Geografici (Minimum Convex Polygons)
# Filtriamo i puntini sulla mappa tenendo solo quelli delle specie condivise
spatial_data <- cleaned_occurrences %>%
  filter(formatted_species %in% shared_species)

# Convertiamo le coordinate in un oggetto "Spaziale" ufficiale per R
sf_points <- st_as_sf(spatial_data, 
                      coords = c("decimalLongitude", "decimalLatitude"), 
                      crs = 4326) # 4326 è il codice standard per le coordinate GPS WGS84

# Raggruppiamo i punti per specie e disegniamo l'areale (Convex Hull)
species_polygons <- sf_points %>%
  group_by(formatted_species) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_convex_hull()

# 8. Visualizzazione Finale!
par(mfrow=c(1,2)) # Dividiamo lo schermo in due
# A sinistra l'albero
plot(pruned_tree, main="Albero Filogenetico (Pruned)")
# A destra i poligoni di overlap ecologico
plot(st_geometry(species_polygons), 
     col = sf.colors(length(shared_species), alpha = 0.5), 
     main="Areali Geografici delle Specie (MCP)")

# 9. Calculate Phylogenetic Distances (WP1 & WP2)
# The cophenetic function extracts the pairwise evolutionary distance between species
phylo_dist_matrix <- cophenetic(pruned_tree)

cat("\n--- PHYLOGENETIC DISTANCES ---\n")
print(phylo_dist_matrix)

# 10. Calculate Geographic Range Areas (WP2 & WP3)
# st_area calculates the area of the polygons in square meters
# We divide by 1,000,000 to convert the result into square kilometers
species_polygons$area_km2 <- as.numeric(st_area(species_polygons)) / 1e6

cat("\n--- GEOGRAPHIC RANGE AREAS (km2) ---\n")
# Print only the species name and its area, dropping the complex geometry column for readability
print(st_drop_geometry(species_polygons)[, c("formatted_species", "area_km2")])

# 11. Optional: Save the results to your computer for the proposal
# write.csv(phylo_dist_matrix, "Phylogenetic_Distances.csv")
# write.csv(st_drop_geometry(species_polygons), "Geographic_Areas.csv")
