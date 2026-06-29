# =============================================================
# 01_build_concatenated_tree.R
# Concatena i 3 allineamenti genici, costruisce l'albero ML/Bayesiano
# e lo time-calibra con un punto di calibrazione fossile/secondario
# =============================================================

# install.packages(c("ape", "phangorn", "ips"))
library(ape)
library(phangorn)

# --- 1. Leggi i tre allineamenti (dopo MAFFT) ----------------------
coi  <- read.phyDat("conus_COI_aligned.fasta", format = "fasta")
s16  <- read.phyDat("conus_16S_aligned.fasta", format = "fasta")
s12  <- read.phyDat("conus_12S_aligned.fasta", format = "fasta")

# --- 2. Verifica quali specie sono presenti in TUTTI i geni --------
# (la tabella supplementare dice che 19 specie hanno solo 1/3 geni:
#  queste vanno escluse dal dataset concatenato, o tenute con gap)

common_species <- Reduce(intersect, list(names(coi), names(s16), names(s12)))
cat("Specie presenti in tutti i 3 geni:", length(common_species), "\n")

coi  <- subset(coi, common_species)
s16  <- subset(s16, common_species)
s12  <- subset(s12, common_species)

# --- 3. Concatena i tre alignment in un'unica matrice --------------
# phangorn non ha una funzione diretta di concatenazione "facile";
# in alternativa si puo' usare il pacchetto 'ips' (mafft + concatenate)
# o costruire la matrice manualmente:

concat_alignment <- cbind(as.character(coi), as.character(s16), as.character(s12))
concat_phydat <- phyDat(concat_alignment, type = "DNA")

write.phyDat(concat_phydat, "conus_concatenated.fasta", format = "fasta")

# --- 4. Costruisci l'albero ML (modello GTR+G, come nel paper originale) ---
dm <- dist.ml(concat_phydat)
tree_nj <- NJ(dm)  # albero iniziale (NJ) come starting tree per ML

fit_start <- pml(tree_nj, data = concat_phydat)
fit_gtr <- update(fit_start, k = 4, inv = 0.2)
fit_gtr <- optim.pml(fit_gtr, model = "GTR", optGamma = TRUE,
                     optInv = TRUE, rearrangement = "stochastic")

tree_ml <- fit_gtr$tree
write.tree(tree_ml, "conus_ML_unrooted.nwk")

cat("\nAlbero ML salvato in conus_ML_unrooted.nwk\n")
cat("Log-likelihood finale:", fit_gtr$logLik, "\n")

# --- 5. Time-calibrazione -------------------------------------------
# Per calibrare in tempo assoluto servono punti di calibrazione
# (fossili o secondari) presi dalla letteratura sui Conidae.
# Due stime indipendenti e coerenti tra loro, da studi con
# calibrazione fossile vera (non solo orologio molecolare puro):
#
#   1. Crown-group CONUS (sensu stretto, dentro il 4-genus split):
#      ~25.7 Ma, 95% HPD: 22.35-28.02 Ma
#      Fonte: Chase et al. 2022, MBE - BEAST2 con calibrazioni fossili,
#      72 taxa, 86 housekeeping genes
#
#   2. Split tra i 4 generi principali (Profundiconus / Californiconus
#      + Lilliconus + Pseudolilliconus / Conasprella / Conus):
#      56-30 Ma (root piu' profondo del nostro albero, se include
#      outgroup di altri generi)
#      Fonte: Uribe, Puillandre & Zardoya 2017 "Beyond Conus", MPE -
#      mitogenomi completi
#
# NOTA IMPORTANTE: il nostro albero (Puillandre et al. 2014) include
# SOLO il genere Conus (sensu lato, con outgroup Conasprella/etc. nelle
# ultime righe della Supplementary Data 2). Quindi il punto di
# calibrazione corretto da usare sul NOSTRO albero e' il #1 (crown
# Conus, ~25.7 Ma), applicato al nodo piu' profondo che separa le
# 4 lineages maggiori (Large Major Clade, Small Major Clade,
# C. californicus, e la quarta lineage identificata da Puillandre).
# Il punto #2 e' piu' profondo (livello di genere) e si applica solo
# se nell'albero sono inclusi anche i 16 outgroup elencati nella
# Supplementary Data 2.

# Identifica il nodo corretto (radice del crown Conus, o radice
# assoluta se sono stati mantenuti gli outgroup):
# plot(tree_ml, cex = 0.4)  # ispeziona visivamente per trovare il nodo
# nodelabels()              # mostra i numeri dei nodi interni

calib <- makeChronosCalib(
  tree_ml,
  node = Ntip(tree_ml) + 1,  # nodo radice (root) - verificare sia il nodo giusto
  age.min = 22.35,
  age.max = 28.02
)

tree_dated <- chronos(tree_ml, calibration = calib, lambda = 1)
write.tree(tree_dated, "conus_timetree.nwk")

cat("\nAlbero time-calibrato salvato in conus_timetree.nwk\n")
cat("Calibrazione applicata: crown Conus 22.35-28.02 Ma (Chase et al. 2022)\n")
cat("Eta' radice nell'albero calibrato:", max(branching.times(tree_dated)), "Ma\n")

# CONTROLLO DI QUALITA' RACCOMANDATO:
# Se nell'albero sono presenti anche gli outgroup (Conasprella,
# Profundiconus, ecc. - controllare gli ultimi nomi della lista in
# Supplementary Data 2), si puo' aggiungere un SECONDO punto di
# calibrazione piu' profondo per il nodo che separa quei generi:
#
# calib2 <- makeChronosCalib(tree_ml, node = <nodo_split_generi>,
#                             age.min = 30, age.max = 56)
# calib_multi <- rbind(calib, calib2)
# tree_dated <- chronos(tree_ml, calibration = calib_multi, lambda = 1)
#
# Usare DUE calibrazioni indipendenti (una piu' superficiale, una piu'
# profonda) riduce l'incertezza e rende la datazione piu' robusta -
# fortemente raccomandato se l'albero include gli outgroup.
