#!/usr/bin/env Rscript

# This script analyzes the diversification history of a specific clade using the relative phylogenetic tree.  
# It generates a Lineage-Through-Time (LTT) plot, computes the Gamma statistic, and calculates the DeltaR statistic.
# The aim of it is to detect signatures of diversification slowdowns or accelerations over time. 
# If the provided tree is heavily undersampled (e.g., 14 out of ~140 known species are present), the script manually performs a Monte Carlo Constant Rates 
# (MCCR) test to correct the Gamma statistic for incomplete taxon sampling, preventing false positive signals of slowdown.
# It requires a specific clade time-calibrated phylogenetic tree in Newick format.


library(ape)

# Import the time-calibrated tree in Newick format (change with the correct name).
tree <- read.tree("NEWICK_TREE")
# Print the total number of tips (species) currently present in the tree.
cat("Number of tips in the tree:", Ntip(tree), "\n")

# Verify if the tree is ultrametric (all tips end at time 0, i.e., the present).
cat("Is the tree ultrametric?", is.ultrametric(tree), "\n")
# Extract branching times and find the maximum to determine the root age (in million years).
cat("Root age (in Myrs):", round(max(branching.times(tree)), 2), "\n")

# Set the plot area to display a single graph.
par(mfrow = c(1, 1))
# Generate the LTT plot with the y-axis on a logarithmic scale.
ltt.plot(tree,
         log = "y",
         main = "LTT plot")

# Extract the exact coordinate points (time and number of lineages) from the plot.
coords <- ltt.plot.coords(tree)
# Save the total number of tips to set the expectation.
n_tips_check <- Ntip(tree)
# Define the x-axis range for the expected constant rate line (from root to present).
expected_line_x <- range(coords[, "time"])
# Define the y-axis range for the expected constant rate line (from 2 lineages to total tips).
expected_line_y <- c(2, n_tips_check)
# Draw a dashed red line representing the expected accumulation under a constant rate.
lines(expected_line_x, expected_line_y, col = "red", lty = 2, lwd = 1.5)
# Add a legend to differentiate the observed data from the theoretical expectation.
legend("topleft", 
       legend = c("Observed", "Constant rate (Expected)"),
       col = c("black", "red"), 
       lty = c(1, 2),
       bty = "n",
       cex = 0.8)

# Calculate the empirical gamma statistic for the tree.
gamma_stat <- gammaStat(tree)
# Print the observed gamma value.
cat("\nObserved Gamma statistic:", round(gamma_stat, 3), "\n")
# Print a reminder on how to interpret the sign of the statistic.
cat("(negative = slowdown, positive = recent acceleration)\n")

# Calculate a naive two-tailed p-value based on a standard normal distribution.
p_naive <- 2 * (1 - pnorm(abs(gamma_stat)))
# Print the naive p-value. It does not account for missing taxa.
cat("Naive p-value (NOT corrected for missing taxa):", round(p_naive, 4), "\n")

# Define the number of species currently in our tree.
n_taxa_tree <- Ntip(tree)
# Define the estimated true total number of species for the clade.
# Change with the correct number. For example, you can check on WoRMS database for marine species.
n_known_taxa <- 140

# Calculate the fraction of the clade that is actually sampled in our tree.
sampling_fraction <- n_taxa_tree / n_known_taxa
# Print the calculated sampling fraction. 
cat("\nSampling fraction:", round(sampling_fraction, 3), "\n")

# Define a function to perform the Monte Carlo Constant Rates (MCCR) test.
# Using 'replicate' replaces the slow 'for' loop and makes the code cleaner.
mccr_manual <- function(n_real, n_sampled, nsim = 1000, seed = 1) {
  # Set the random seed to ensure reproducibility of the simulation.
  set.seed(seed)
  # Repeat the enclosed expression 'nsim' times and store results in a vector.
  replicate(nsim, {
    # Simulate a full tree with 'n_real' species under a pure birth model.
    full_tree <- rphylo(n = n_real, birth = 1, death = 0)
    # Randomly select a subset of tips matching our empirical sample size ('n_sampled').
    pruned_tips <- sample(full_tree$tip.label, n_sampled)
    # Prune the simulated tree to keep only the randomly sampled tips.
    pruned_tree <- keep.tip(full_tree, pruned_tips)
    # Calculate and return the gamma statistic for this pruned simulated tree.
    gammaStat(pruned_tree)
  })
}

# Execute the MCCR function to generate a null distribution of gamma values.
gamma_null_dist <- mccr_manual(n_known_taxa, n_taxa_tree, nsim = 1000)

# Calculate the corrected p-value: proportion of simulated gammas as or more extreme than the observed one.
p_mccr <- mean(abs(gamma_null_dist) >= abs(gamma_stat))
# Calculate the critical threshold (bottom 2.5% tail) of the simulated null distribution.
critical_value <- quantile(gamma_null_dist, 0.025)

# Print the mean of the simulated null distribution.
cat("  Mean:", round(mean(gamma_null_dist), 3), "\n")
# Print the critical threshold boundary.
cat("  Critical threshold (2.5th percentile):", round(critical_value, 3), "\n")
# Reprint the observed empirical gamma for direct visual comparison.
cat("Observed Gamma:", round(gamma_stat, 3), "\n")
# Print the final, reliable, MCCR-corrected two-tailed p-value.
cat("MCCR p-value (two-tailed):", round(p_mccr, 4), "\n")


# Define a function to calculate an alternative diversification metric (DeltaR).
deltaR_simple <- function(tree) {
  # Extract the LTT plot coordinates again.
  coords <- ltt.plot.coords(tree)
  # Calculate the total evolutionary time spanned by the tree (root to tip).
  total_time <- max(coords[, "time"]) - min(coords[, "time"])
  # Determine the temporal midpoint of the tree.
  midpoint <- min(coords[, "time"]) + total_time / 2
  
  # Subset coordinates representing the first half of the tree's history.
  early <- coords[coords[, "time"] <= midpoint, ]
  # Subset coordinates representing the second half of the tree's history.
  late  <- coords[coords[, "time"] >  midpoint, ]
  
  # Calculate diversification rate for the early half (log(lineages ratio) / time span).
  rate_early <- log(max(early[, "N"]) / min(early[, "N"])) / (max(early[,"time"]) - min(early[,"time"]))
  # Calculate diversification rate for the late half (log(lineages ratio) / time span).
  rate_late  <- log(max(late[, "N"])  / min(late[, "N"]))  / (max(late[,"time"])  - min(late[,"time"]))
  
  # Return a list containing early rate, late rate, and the difference (DeltaR).
  list(rate_early = rate_early, rate_late = rate_late, deltaR = rate_late - rate_early)
}

# Run the DeltaR calculation on our empirical tree.
deltaR_res <- deltaR_simple(tree)

# Print the rate calculated for the first half of the tree.
cat("Rate first half:", round(deltaR_res$rate_early, 4), "\n")
# Print the rate calculated for the second half of the tree.
cat("Rate second half:", round(deltaR_res$rate_late, 4), "\n")
# Print the final DeltaR difference (negative implies a slowdown in the second half).
cat("DeltaR:", round(deltaR_res$deltaR, 4), " (negative = slowdown)\n")

# Create a single-row dataframe consolidating all the calculated metrics for the clade.
summary_row <- data.frame(
  clade = "CLADE_NAME", # Change with the correct name.
  n_tips = n_taxa_tree,
  n_known_species = n_known_taxa,
  sampling_fraction = sampling_fraction,
  gamma = gamma_stat,
  p_naive = p_naive,
  p_mccr = p_mccr,
  deltaR = deltaR_res$deltaR
)

# Export the summary row to a CSV file (change the file name as desired).
write.csv(summary_row, "DIVERISFICATION_SUMMARY.CSV", row.names = FALSE)
print(summary_row)
