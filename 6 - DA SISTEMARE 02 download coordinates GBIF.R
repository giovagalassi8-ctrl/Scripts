# =============================================================
# 02_download_coordinates_Carditidae.R
# Scarica le coordinate GBIF per le 14 specie di Carditidae,
# usando i nomi accettati risolti al passo 1
# =============================================================

library(rgbif)
library(dplyr)

# --- 1. Leggi la tabella di lookup tassonomico dal passo 1 ----------
lookup <- read.csv("Carditidae_taxonomy_lookup_worms.csv", stringsAsFactors = FALSE)

# Usa valid_name se disponibile, altrimenti il nome originale
lookup$name_to_query <- ifelse(!is.na(lookup$valid_name) & lookup$valid_name != "",
                               lookup$valid_name, lookup$queried_name)

cat("Nomi che verranno usati per la query GBIF:\n")
print(lookup$name_to_query)

# --- 2. Scarica le occorrenze GBIF per ciascuna specie --------------
download_gbif <- function(species_name) {
  cat("Scaricando:", species_name, "...\n")
  res <- tryCatch(
    occ_search(scientificName = species_name, hasCoordinate = TRUE, limit = 2000),
    error = function(e) {
      cat("  Errore per", species_name, ":", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(res) || is.null(res$data) || nrow(res$data) == 0) {
    cat("  Nessuna occorrenza trovata per", species_name, "\n")
    return(NULL)
  }
  res$data
}

all_occurrences <- lapply(lookup$name_to_query, download_gbif)
all_occurrences <- all_occurrences[!sapply(all_occurrences, is.null)]

# --- 3. Combina tutto in un unico dataframe -------------------------
# Usiamo solo le colonne essenziali per evitare conflitti di tipo
# tra le risposte di specie diverse
essential_cols <- c("scientificName", "species", "decimalLatitude",
                    "decimalLongitude", "depth", "country", "year")

combined <- lapply(all_occurrences, function(df) {
  cols_present <- intersect(essential_cols, colnames(df))
  df[, cols_present, drop = FALSE]
})
combined_df <- bind_rows(combined)

cat("\nTotale occorrenze scaricate:", nrow(combined_df), "\n")
cat("Specie con almeno un'occorrenza:", length(unique(combined_df$species)), "/", nrow(lookup), "\n")

# --- 4. Controllo qualita': occorrenze per specie --------------------
occ_per_species <- combined_df %>% count(species, name = "n_occurrences") %>% arrange(n_occurrences)
print(occ_per_species)

MIN_OCC <- 5
sp_pochi_dati <- occ_per_species$species[occ_per_species$n_occurrences < MIN_OCC]
if (length(sp_pochi_dati) > 0) {
  cat("\nATTENZIONE - specie con meno di", MIN_OCC, "occorrenze (overlap poco affidabile):\n")
  print(sp_pochi_dati)
}

# Specie completamente assenti da GBIF
missing_species <- setdiff(lookup$name_to_query, unique(combined_df$species))
if (length(missing_species) > 0) {
  cat("\nATTENZIONE - specie SENZA NESSUNA occorrenza GBIF:\n")
  print(missing_species)
}

# --- 5. Salva il file finale -----------------------------------------
write.csv(combined_df, "Carditidae_coordinates_resolved.csv", row.names = FALSE)
cat("\nSalvato in Carditidae_coordinates_resolved.csv\n")
cat("Questo e' il file da usare nello script di overlap geografico.\n")
