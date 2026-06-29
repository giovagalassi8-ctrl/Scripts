# =============================================================
# 01_resolve_taxonomy_Carditidae.R
# Risoluzione tassonomica dei 14 tip di Carditidae (WoRMS + GBIF fallback)
# =============================================================

# install.packages(c("worrms", "rgbif", "stringr"))
library(ape)
library(worrms)
library(rgbif)
library(stringr)

# --- 1. Leggi l'albero potato (solo Carditidae) --------------------
tree <- read.tree("Carditidae.nwk")
tip_names <- tree$tip.label
tip_species <- str_replace_all(tip_names, "_", " ")

cat("Numero di tip nell'albero:", length(tip_species), "\n")
print(tip_species)

# --- 2. Funzione di risoluzione via WoRMS -----------------------
resolve_worms <- function(sp_name) {
  res <- tryCatch(
    wm_records_taxamatch(name = sp_name, marine_only = TRUE),
    error = function(e) NULL
  )
  if (is.null(res) || length(res) == 0 || nrow(res[[1]]) == 0) {
    return(data.frame(
      queried_name = sp_name, status = "NOT_FOUND",
      valid_name = NA, valid_AphiaID = NA
    ))
  }
  rec <- res[[1]][1, ]
  data.frame(
    queried_name = sp_name,
    status = rec$status,
    valid_name = rec$valid_name,
    valid_AphiaID = rec$valid_AphiaID
  )
}

# --- 3. Applica a tutti i tip -----------------------------------
worms_results <- do.call(rbind, lapply(tip_species, resolve_worms))
print(worms_results)

# --- 4. Per i NOT_FOUND, fallback su GBIF -----------------------
not_found <- worms_results$queried_name[worms_results$status == "NOT_FOUND"]

if (length(not_found) > 0) {
  cat("\nSpecie non trovate su WoRMS, provo fallback GBIF:\n")
  print(not_found)
  gbif_fallback <- lapply(not_found, function(nm) {
    res <- name_backbone(name = nm, rank = "species")
    data.frame(
      queried_name = nm,
      gbif_matchType = ifelse(is.null(res$matchType), NA, res$matchType),
      gbif_status = ifelse(is.null(res$status), NA, res$status),
      gbif_accepted_name = ifelse(is.null(res$species), NA, res$species),
      gbif_usageKey = ifelse(is.null(res$usageKey), NA, res$usageKey)
    )
  })
  gbif_fallback <- do.call(rbind, gbif_fallback)
  print(gbif_fallback)
  write.csv(gbif_fallback, "Carditidae_taxonomy_gbif_fallback.csv", row.names = FALSE)
}

# --- 5. Salva la tabella di lookup finale ------------------------
write.csv(worms_results, "Carditidae_taxonomy_lookup_worms.csv", row.names = FALSE)

cat("\n--- RIEPILOGO ---\n")
cat("Trovati su WoRMS:", sum(worms_results$status != "NOT_FOUND"), "/", length(tip_species), "\n")
cat("Da controllare manualmente (NOT_FOUND):", sum(worms_results$status == "NOT_FOUND"), "\n")

# NOTA: con solo 14 specie, il controllo manuale della tabella
# finale richiede pochi minuti - fallo comunque prima di scaricare
# le coordinate, controllando che valid_name sia sensato per ognuna.
