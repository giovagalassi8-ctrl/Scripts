# =============================================================
# 02b_calibrate_no_outgroups.R
# Time-calibrazione dell'albero IQ-TREE costruito SENZA outgroup.
# La radice dell'albero E' GIA' il crown-group di Conus - un solo
# punto di calibrazione, robusto, niente problemi di monofilia.
# =============================================================

library(ape)

# --- 1. Leggi l'albero ML (costruito su conus_concatenated.fasta SENZA
# outgroup - il nome del file di output mantiene "conus_concatenated"
# perche' IQ-TREE usa il nome del file di input come prefisso, anche se
# il contenuto del file di input era gia' stato pulito dagli outgroup)
tree_ml <- read.tree("conus_concatenated.fasta.treefile")
cat("Numero di tip nell'albero:", Ntip(tree_ml), "\n")
stopifnot(Ntip(tree_ml) == 247)  # controllo di sicurezza: deve essere
# il file SENZA outgroup (247 tip),
# non quello vecchio CON outgroup (260)

# --- 2. Controllo: l'albero deve essere binario (anche se unrooted, --
# per chronos() la trifurcazione alla radice non e' un problema in se',
# ma e' piu' sicuro verificare la struttura prima di procedere)
cat("Albero binario?", is.binary(tree_ml), "\n")
if (!is.binary(tree_ml)) {
  tree_ml <- multi2di(tree_ml)
  cat("Politomie risolte con multi2di()\n")
}

# --- 3. Calibrazione: crown-group Conus ~25.7 Ma --------------------
# Fonte: Chase et al. 2022, MBE - BEAST2 con calibrazioni fossili
# 95% HPD: 22.35 - 28.02 Ma
# Applicata sulla radice assoluta, che ora E' il crown-Conus
# (non avendo piu' outgroup nell'albero)

root_node <- Ntip(tree_ml) + 1

calib <- makeChronosCalib(
  tree_ml,
  node = root_node,
  age.min = 22.35,
  age.max = 28.02
)
print(calib)

# --- 4. Esegui la time-calibrazione (penalized likelihood) ----------
tree_dated <- chronos(tree_ml, calibration = calib, lambda = 1)

cat("\nAlbero time-calibrato creato.\n")
cat("Eta' radice:", round(max(branching.times(tree_dated)), 2), "Ma\n")

# --- 5. Salva il risultato -------------------------------------------
write.tree(tree_dated, "conus_timetree_final.nwk")
cat("Salvato in conus_timetree_final.nwk - QUESTO e' il file da usare\n")
cat("per gli script successivi (LTT, gamma, overlap geografico)\n")

cat("\n========== RIEPILOGO ==========\n")
cat("Tip nel time-tree finale (crown Conus): ", Ntip(tree_dated), "\n")
cat("Eta' crown-group:                       ", round(max(branching.times(tree_dated)), 2), "Ma\n")
cat("================================\n")

# --- 6. Controllo visivo finale --------------------------------------
plot(tree_dated, cex = 0.4, no.margin = TRUE)
axisPhylo()

# NOTE:
# - Questo approccio usa un solo punto di calibrazione (non la doppia
#   calibrazione che avevamo tentato con gli outgroup) - e' meno
#   "ridondante" ma molto piu' robusto, perche' evita completamente
#   il problema di monofilia incerta riscontrato con gli outgroup.
# - Se in futuro si vuole rafforzare la calibrazione con un secondo
#   punto INTERNO (es. l'eta' di uno split noto tra sottogeneri di
#   Conus, se disponibile in letteratura), si puo' aggiungere come
#   seconda riga a calib, identificando il nodo con getMRCA() su un
#   sottoinsieme di specie note per quel clade specifico.
