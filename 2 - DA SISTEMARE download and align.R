# =============================================================
# 00_download_and_align.R
# Scarica le sequenze GenBank (COI, 16S, 12S) per Conus usando
# gli accession number di Supplementary Data 2 (Puillandre et al. 2014)
# e prepara l'allineamento per ciascun gene
# =============================================================

# install.packages(c("rentrez", "ape", "msa", "readxl"))
# Bioconductor: BiocManager::install("msa")
library(rentrez)
library(ape)
library(readxl)

# --- 1. Leggi la tabella Supplementary Data 2 ---------------------
# Apri il file mmc2.docx convertendolo prima in .csv/.xlsx a mano
# (e' un Word document - copia la tabella in Excel/Google Sheets e
#  esporta come .csv con colonne: species, voucher, prey, province,
#  accession_COI, accession_16S, accession_12S)

sp_table <- read.csv("Conus.csv", stringsAsFactors = FALSE)

cat("Numero di specie nella tabella:", nrow(sp_table), "\n")
# Dovrebbero essere 345 - 19 (solo 1 gene) - 16 (outgroup) = ~310 utilizzabili
# per il dataset concatenato completo

# --- 2. Funzione per scaricare sequenze da GenBank -----------------
# IMPORTANTE: la funzione ritorna una lista NOMINATA per species,
# non per accession - cosi' il nome resta sempre agganciato alla
# sequenza corretta anche se alcuni accession sono mancanti/vuoti
# o se il download fallisce per qualcuno di essi (niente disallineamento
# di indici tra sequenze e nomi specie).

download_seqs <- function(accessions, species_names, gene_name) {
  valid <- !is.na(accessions) & accessions != ""
  cat("Scaricando", sum(valid), "sequenze per", gene_name,
      "(", sum(!valid), "specie senza accession per questo gene )\n")
  
  seqs <- setNames(vector("list", length(accessions)), species_names)
  
  for (i in seq_along(accessions)) {
    if (!valid[i]) next  # specie senza accession per questo gene: resta NULL
    acc <- accessions[i]
    result <- tryCatch({
      entrez_fetch(db = "nuccore", id = acc, rettype = "fasta")
    }, error = function(e) {
      cat("  Errore per accession", acc, "(", species_names[i], "):",
          conditionMessage(e), "\n")
      NULL
    })
    seqs[[species_names[i]]] <- result
    Sys.sleep(0.34)  # rispetta i rate limit NCBI (max ~3 richieste/secondo)
  }
  return(seqs)
}

# --- 3. Scarica per ciascun gene -----------------------------------
coi_seqs <- download_seqs(sp_table$accession_COI, sp_table$species, "COI")
s16_seqs <- download_seqs(sp_table$accession_16S, sp_table$species, "16S")
s12_seqs <- download_seqs(sp_table$accession_12S, sp_table$species, "12S")

# --- 4. Salva in formato FASTA -------------------------------------
# Scrive una sequenza per specie SOLO se il download e' andato a buon
# fine (seq_list[[nome]] non e' NULL) - usa il nome della lista come
# header, non la posizione, quindi non puo' disallinearsi
write_fasta <- function(seq_list, filename) {
  fasta_lines <- c()
  n_written <- 0
  for (sp_name in names(seq_list)) {
    if (!is.null(seq_list[[sp_name]])) {
      header <- paste0(">", sp_name)
      seq_body <- sub("^>[^\n]*\n", "", seq_list[[sp_name]])  # rimuove header GenBank originale
      fasta_lines <- c(fasta_lines, header, seq_body)
      n_written <- n_written + 1
    }
  }
  writeLines(fasta_lines, filename)
  cat("Scritte", n_written, "sequenze in", filename, "\n")
}

write_fasta(coi_seqs, "conus_COI_raw.fasta")
write_fasta(s16_seqs, "conus_16S_raw.fasta")
write_fasta(s12_seqs, "conus_12S_raw.fasta")

cat("\nFile FASTA salvati. Prossimo passo: allineamento con MAFFT.\n")
cat("Da terminale (fuori da R):\n")
cat("  mafft --auto conus_COI_raw.fasta > conus_COI_aligned.fasta\n")
cat("  mafft --auto conus_16S_raw.fasta > conus_16S_aligned.fasta\n")
cat("  mafft --auto conus_12S_raw.fasta > conus_12S_aligned.fasta\n")

# NOTE:
# - rentrez fa query a NCBI Entrez via API pubblica e gratuita,
#   ma e' buona norma registrare una email/API key (set_entrez_key())
#   per rate limit piu' alti se le sequenze sono centinaia.
# - Se alcuni accession non si scaricano (sequenza ritirata/cambiata
#   da quando il paper e' stato pubblicato nel 2014), e' normale -
#   segnare quali specie risultano mancanti e procedere con le altre.
