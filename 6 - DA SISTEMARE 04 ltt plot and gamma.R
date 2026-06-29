# =============================================================
# 04_ltt_and_gamma_Carditidae.R
# LTT plot + gamma statistic (Pybus & Harvey) + correzione MCCR
# per il clade Carditidae (14 tip)
# =============================================================

library(ape)

# --- 1. Leggi l'albero -----------------------------------------------
tree <- read.tree("Carditidae.nwk")
cat("Numero di tip nell'albero:", Ntip(tree), "\n")

# IMPORTANTE: questo albero NON e' ancora time-calibrato (i branch
# length vengono direttamente da timetree.org, che dovrebbe già
# fornirli in milioni di anni - verifichiamolo)
cat("Albero ultrametrico?", is.ultrametric(tree), "\n")
cat("Eta' radice (se in milioni di anni):", round(max(branching.times(tree)), 2), "\n")

# Se is.ultrametric() e' FALSE, l'albero ha qualche inconsistenza nei
# branch length (frequente quando si pota un albero piu' grande) -
# in tal caso usare force.ultrametric() di phytools come correzione
# leggera, oppure ri-verificare l'albero originale prima di procedere.

# --- 2. LTT plot ------------------------------------------------------
par(mfrow = c(1, 1))
ltt.plot(tree, log = "y", main = "LTT plot - Carditidae")

coords <- ltt.plot.coords(tree)
n_tips_check <- Ntip(tree)
expected_line_x <- range(coords[, "time"])
expected_line_y <- c(2, n_tips_check)
lines(expected_line_x, expected_line_y, col = "red", lty = 2, lwd = 1.5)
legend("topleft", legend = c("Osservato", "Tasso costante (atteso)"),
       col = c("black", "red"), lty = c(1, 2), bty = "n", cex = 0.8)

# --- 3. Gamma statistic (Pybus & Harvey 1995) -------------------------
gamma_stat <- gammaStat(tree)
cat("\nGamma statistic osservato:", round(gamma_stat, 3), "\n")
cat("(negativo = slowdown, positivo = accelerazione recente)\n")

p_naive <- 2 * (1 - pnorm(abs(gamma_stat)))
cat("p-value naive (NON corretto per missing taxa):", round(p_naive, 4), "\n")

# --- 4. Correzione MCCR per sotto-campionamento ------------------------
# NOTA SUL DENOMINATORE: la famiglia Carditidae conta circa 140
# specie viventi note (WoRMS) - il nostro albero ne campiona solo 14,
# quindi la sampling fraction e' molto bassa (~10%). Questo era
# previsto fin dall'inizio (vedi discussione iniziale sulla
# fattibilita') ed e' uno dei limiti piu' importanti da dichiarare
# esplicitamente per questo specifico clade.

n_taxa_albero <- Ntip(tree)             # 14
n_taxa_noti <- 140                      # Carditidae, sensu WoRMS

sampling_fraction <- n_taxa_albero / n_taxa_noti
cat("\nFrazione campionata (14/140):", round(sampling_fraction, 3), "\n")

# Implementazione diretta del MCCR test (vedi nota nello script di
# Conus sul perche' evitiamo geiger::mccr - conflitti di namespace
# riscontrati in precedenza)
mccr_manual <- function(n_real, n_sampled, nsim = 1000, seed = 1) {
  set.seed(seed)
  gamma_null <- numeric(nsim)
  for (i in 1:nsim) {
    full_tree <- rphylo(n = n_real, birth = 1, death = 0)
    pruned_tips <- sample(full_tree$tip.label, n_sampled)
    pruned_tree <- keep.tip(full_tree, pruned_tips)
    gamma_null[i] <- gammaStat(pruned_tree)
  }
  gamma_null
}

cat("\nEseguo 1000 simulazioni MCCR...\n")
gamma_null_dist <- mccr_manual(n_taxa_noti, n_taxa_albero, nsim = 1000)

p_mccr <- mean(abs(gamma_null_dist) >= abs(gamma_stat))
critical_value <- quantile(gamma_null_dist, 0.025)

cat("\nDistribuzione nulla simulata (sotto sotto-campionamento casuale):\n")
cat("  Media:", round(mean(gamma_null_dist), 3), "\n")
cat("  Soglia critica (2.5° percentile):", round(critical_value, 3), "\n")
cat("Gamma osservato:", round(gamma_stat, 3), "\n")
cat("p-value MCCR (due code):", round(p_mccr, 4), "\n")

# --- 5. DeltaR statistic ------------------------------------------------
deltaR_simple <- function(tree) {
  coords <- ltt.plot.coords(tree)
  total_time <- max(coords[, "time"]) - min(coords[, "time"])
  midpoint <- min(coords[, "time"]) + total_time / 2
  early <- coords[coords[, "time"] <= midpoint, ]
  late  <- coords[coords[, "time"] >  midpoint, ]
  rate_early <- log(max(early[, "N"]) / min(early[, "N"])) / (max(early[,"time"]) - min(early[,"time"]))
  rate_late  <- log(max(late[, "N"])  / min(late[, "N"]))  / (max(late[,"time"])  - min(late[,"time"]))
  list(rate_early = rate_early, rate_late = rate_late, deltaR = rate_late - rate_early)
}

deltaR_res <- deltaR_simple(tree)
cat("\n--- DeltaR ---\n")
cat("Tasso prima meta':", round(deltaR_res$rate_early, 4), "\n")
cat("Tasso seconda meta':", round(deltaR_res$rate_late, 4), "\n")
cat("DeltaR:", round(deltaR_res$deltaR, 4), " (negativo = slowdown)\n")

# --- 6. Salva il riepilogo -----------------------------------------------
summary_row <- data.frame(
  clade = "Carditidae",
  n_tips = n_taxa_albero,
  n_known_species = n_taxa_noti,
  sampling_fraction = sampling_fraction,
  gamma = gamma_stat,
  p_naive = p_naive,
  p_mccr = p_mccr,
  deltaR = deltaR_res$deltaR
)
write.csv(summary_row, "diversification_summary_Carditidae.csv", row.names = FALSE)
print(summary_row)

# NOTA: con solo 14 tip su 140 specie note (~10% sampling), questo
# risultato va interpretato con grande cautela - utile come punto
# nel confronto aggregato multi-clade, non come stima robusta da sola.
