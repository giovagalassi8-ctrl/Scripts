# =============================================================
# 03_geographic_overlap_Carditidae.R
# Calcola l'overlap geografico tra coppie di specie di Carditidae
# Metodo: convex hull (semplice) + Schoener's D via kernel density
# =============================================================

# install.packages(c("sf", "dplyr", "ggplot2"))
library(sf)
library(dplyr)
library(ggplot2)

# --- 1. Carica i dati GBIF (gia' filtrati per famiglia/clade) --
# Assumiamo un CSV con almeno: species, decimalLatitude, decimalLongitude
# e che i nomi siano GIA' stati armonizzati con lo step di tassonomia (punto 3)

occ <- read.csv("Carditidae_coordinates_resolved.csv", stringsAsFactors = FALSE)

# pulizia base: rimuovi righe senza coordinate o coordinate a (0,0)
occ <- occ %>%
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>%
  filter(!(decimalLatitude == 0 & decimalLongitude == 0))

# quante occorrenze per specie? (utile per decidere soglie minime)
occ_per_species <- occ %>% count(species, name = "n_occurrences") %>% arrange(n_occurrences)
print(occ_per_species)

# soglia minima ragionevole per costruire un range attendibile:
# sotto 5 punti il convex hull e' instabile -> segnaliamo ma non escludiamo
# automaticamente (decisione da prendere insieme a Giovanni)
MIN_OCC <- 5
sp_pochi_dati <- occ_per_species$species[occ_per_species$n_occurrences < MIN_OCC]
if (length(sp_pochi_dati) > 0) {
  cat("ATTENZIONE - specie con meno di", MIN_OCC, "occorrenze (overlap poco affidabile):\n")
  print(sp_pochi_dati)
}

# --- 2. Costruisci un poligono (convex hull) per ogni specie -----
# Converti in oggetto sf
occ_sf <- st_as_sf(occ, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

build_hull <- function(sp_name, data_sf) {
  pts <- data_sf %>% filter(species == sp_name)
  if (nrow(pts) < 3) {
    return(NULL)  # serve almeno un triangolo per un poligono
  }
  hull <- st_convex_hull(st_union(pts))
  return(hull)
}

species_list <- unique(occ$species)
hulls <- lapply(species_list, build_hull, data_sf = occ_sf)
names(hulls) <- species_list
hulls <- hulls[!sapply(hulls, is.null)]  # rimuovi specie senza hull valido

cat("\nHull costruiti per", length(hulls), "specie su", length(species_list), "\n")

# --- 3. Calcola l'overlap pairwise (area condivisa / area unione) ----
# Questo e' un indice semplice di overlap (Jaccard-like su poligoni)
# In alternativa a un kernel-density Schoener's D, e' molto piu' robusto
# quando i dati sono pochi - consigliato come prima passata

pairwise_overlap <- function(hull_list) {
  sp_names <- names(hull_list)
  n <- length(sp_names)
  out <- data.frame()
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      h1 <- hull_list[[i]]
      h2 <- hull_list[[j]]
      inter_area <- tryCatch(
        st_area(st_intersection(h1, h2)),
        error = function(e) 0
      )
      union_area <- tryCatch(
        st_area(st_union(h1, h2)),
        error = function(e) NA
      )
      if (length(inter_area) == 0) inter_area <- 0
      overlap_idx <- as.numeric(inter_area) / as.numeric(union_area)
      out <- rbind(out, data.frame(
        sp1 = sp_names[i], sp2 = sp_names[j],
        overlap_index = overlap_idx
      ))
    }
  }
  return(out)
}

overlap_table <- pairwise_overlap(hulls)
print(overlap_table)

write.csv(overlap_table, "pairwise_overlap_Carditidae.csv", row.names = FALSE)

# --- 4. Sintetizza un singolo valore di overlap per il CLADE ------
# Media di tutti gli overlap pairwise = indice di overlap del clade
# (questo e' il numero che andra' sull'asse X del grafico finale
#  overlap vs slowdown, uno per ciascun clade/famiglia)

clade_overlap_mean <- mean(overlap_table$overlap_index, na.rm = TRUE)
cat("\nOverlap medio del clade Carditidae:", round(clade_overlap_mean, 3), "\n")

# --- 5. Visualizzazione rapida (utile per controllo qualita') ----
all_hulls_sf <- do.call(c, hulls)
plot(st_geometry(all_hulls_sf), col = adjustcolor(1:length(hulls), alpha.f = 0.3),
     border = 1:length(hulls), main = "Range geografici (convex hull) - Carditidae")

# NOTE METODOLOGICHE:
# - Il convex hull sovrastima il range reale (include aree mai
#   occupate dentro il poligono). Per un'analisi piu' raffinata,
#   in futuro si puo' passare ad alpha-hull (package 'rangeBuilder')
#   o kernel density (package 'adehabitatHR', metodo Schoener's D).
# - Per ora, dato che e' un proof-of-concept, il convex hull e'
#   accettabile e molto piu' rapido da implementare su molte specie.
