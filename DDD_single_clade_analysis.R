#!/usr/bin/env Rscript

# This script performs a Diversity-Dependent Diversification (DDD) analysis on one clade at a time.
# It fits and compares three diversification models on the time-calibrated tree of the selected clade:
#   CR   = Constant Rate: speciation and extinction rates fixed over time (null model, 2 parameters)
#   DD_L = the speciation rate DECLINES as the number of lineages approaches the carrying capacity K (3 parameters)
#   DD_M = the extinction rate INCREASES as the number of lineages approaches the carrying capacity K (3 parameters)
# The models are compared using the Akaike Information Criterion (AIC) and the Akaike weights.
# It requires a data folder containing one sub-folder per clade, each one with an ultrametric tree in Newick format (.nwk),
# and the total number of accepted extant species of the clade (taken from online databases or from literature).
# The folder of the selected clade is found automatically by searching the data folder for a sub-folder matching the clade name, whatever phylum it sits under.
# To analyse another clade, change the clade name and its known species count, and run the script again
# (use Run / Ctrl+Enter instead of source(), so that the progress messages appear as they happen).
# Results accumulate: each run updates a combined CSV file holding every clade analysed so far, so there is no need to re-run the previous ones.
# Outputs (saved in the output directory):
#   <clade>_dd_fits.rds        = raw model fits for the selected clade
#   <clade>_dd_comparison.csv  = AIC comparison table for the selected clade
#   ALL_clades_comparison.csv  = combined table of all the clades, updated after every run


library(ape)
library(DDD)
library(dplyr)

# Set the clade to analyse (change accordingly).
# Expected structure of the repository: there are several numbered sub-folders of which one, usually the first, contains the data (in this case, it is called 00_Data).
# The data folder contains other sub-folders, one per clade, containing all the data for each group, including a treefile in Newick format for that group.
# The name selected below must match the clade data folder name exactly.
TARGET_CLADE <- "CLADE_NAME"

# Set the total number of accepted extant species of the clade (taken from online databases, like WoRMS, or from literature; change accordingly).
KNOWN_SPECIES <- 1000

# Set the data folder (change accordingly).
DATA_ROOT <- "DATA_FOLDER_NAME"
# Set the output directory name (change accordingly).
output_dir     <- "OUTPUT_DIRECTORY_NAME"

# Create the output directory if it does not exist.
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Search every sub-folder in the data folder and keep the one whose name matches the clade.
all_dirs <- list.dirs(DATA_ROOT, recursive = TRUE, full.names = TRUE)
clade_folder <- all_dirs[basename(all_dirs) == TARGET_CLADE]

# Stop the process if the selected name does not match the name of any sub-folder.
if (length(clade_folder) == 0)
  stop("No folder named '", TARGET_CLADE, "' found under ", DATA_ROOT)
# Stop the process if there are multiple sub-folders with the same matching name.
if (length(clade_folder) > 1)
  stop("Multiple folders named '", TARGET_CLADE, "' found - resolve the ambiguity.")

# Look for the Newick treefile in the folder of the selected clade.
# Stop the process if no Newick treefile is found.
nwk_files <- list.files(clade_folder, pattern = "\\.nwk$", full.names = TRUE)
if (length(nwk_files) == 0)
  stop("No .nwk file found in ", clade_folder)
# Read the treefile and check if it contains more than one tree (in this case only the first one is used) and if the tree is ultrametric.
tree <- read.tree(nwk_files[1])
if (inherits(tree, "multiPhylo")) tree <- tree[[1]]
if (!is.ultrametric(tree))
  stop("Tree is not ultrametric - DDD requires an ultrametric tree.")

# Create an object containing the number of tips of the tree.
n_tips <- Ntip(tree)
# Calculate the number of missing species comparing the number of tips in the tree with the number of accepted extant species set before.
missnumspec <- max(0, KNOWN_SPECIES - n_tips)
# Create an object containing the branching times of the tree.
brts <- branching.times(tree)
# Print a summary of some tree statistics as screen output.
# In particular, the ratio between missing species and sampled tips is returned: with large ratio values (especially greater than 12)
# the script may take much longer, as the model is more complex to manage.
cat("Tips:", n_tips, "| Known species:", KNOWN_SPECIES, "| Missing:", missnumspec, "| Ratio:", round(missnumspec / n_tips, 1), "\n")

# Warn if the tree has more tips than the number of known species.
if (n_tips > KNOWN_SPECIES)
  warning("Tree has MORE tips than the stated known species count.")
# Check if the ratio is higher than 12 and, in case, warn that the process could require a lot of time and that K could be poorly constrained.
if (missnumspec / n_tips > 12)
  cat("WARNING: missing species greatly outnumber sampled tips.",
      "This fit may be slow and K poorly constrained.\n")

# Derive the starting values of the optimisation from the tree itself.
# A fixed starting guess (e.g. lambda = 0.3 for every clade) can fail on large or fast-diversifying trees: the likelihood at that starting point
# can underflow to -Inf, or the optimiser can wander into a degenerate loop (as happened with Dugesiidae, 215 tips).
# Deriving the starting values from the tree keeps the optimiser in a plausible region for the actual scale of the selected clade,
# without changing the true maximum likelihood: it only changes where the search begins.
# OPTIONALLY: to use the same fixed starting values for every clade, uncomment the two following lines and comment the derived 'initial_lambda' and 'initial_mu' below (change accordingly).
#initial_lambda <- 0.3
#initial_mu     <- 0.1

# Calculate the crown age of the tree (the oldest branching time).
tree_age       <- max(brts)
# Set the initial lambda value as a rough estimate of the net diversification rate: ln(number of tips) / crown age.
initial_lambda <- log(n_tips) / tree_age
# Set the initial mu value as a fraction of the speciation rate (it is only a starting guess).
initial_mu     <- initial_lambda * 0.3

# Print the derived starting values.
cat("Derived starting values: lambda =", round(initial_lambda, 4),
    "| mu =", round(initial_mu, 4), "\n")

# Condition the likelihood on the survival of the clade.
cond      <- 1
# Set the initial K value above the total diversity of the clade (sampled tips + missing species).
initial_K <- (n_tips + missnumspec) * 2

# Fit the Constant Rate (CR) birth-death model, estimating lambda and mu.
# If a fit fails, the error message is printed and the result is set to NULL, so that the script can continue with the other models.
cat("\nFitting CR...\n")
cr_fit <- tryCatch(
  bd_ML(brts = brts, initparsopt = c(initial_lambda, initial_mu),
        idparsopt = c(1, 2), missnumspec = missnumspec, cond = cond),
  error = function(e) { cat("  failed:", e$message, "\n"); NULL })

# Fit the DD_L model (ddmodel = 1: linear diversity-dependence of the speciation rate), estimating lambda, mu and K.
cat("Fitting DD_L (diversity-dependent speciation)...\n")
ddl_fit <- tryCatch(
  dd_ML(brts = brts, initparsopt = c(initial_lambda, initial_mu, initial_K),
        idparsopt = c(1, 2, 3), missnumspec = missnumspec, cond = cond,
        ddmodel = 1),
  error = function(e) { cat("  failed:", e$message, "\n"); NULL })

# Fit the DD_M model (ddmodel = 3: linear diversity-dependence of the extinction rate), estimating lambda, mu and K.
cat("Fitting DD_M (diversity-dependent extinction)...\n")
ddm_fit <- tryCatch(
  dd_ML(brts = brts, initparsopt = c(initial_lambda, initial_mu, initial_K),
        idparsopt = c(1, 2, 3), missnumspec = missnumspec, cond = cond,
        ddmodel = 3),
  error = function(e) { cat("  failed:", e$message, "\n"); NULL })

# Save the raw model fits of the selected clade, together with its main information, in .rds format.
saveRDS(list(clade = TARGET_CLADE, n_tips = n_tips, missnumspec = missnumspec,
             known_species = KNOWN_SPECIES,
             CR = cr_fit, DD_L = ddl_fit, DD_M = ddm_fit),
        file.path(output_dir, paste0(TARGET_CLADE, "_dd_fits.rds")))

# Define a function that builds the row of the comparison table for a model, extracting the estimated parameters, the log-likelihood and the AIC.
# If the model failed to fit, the row is filled with NA values and flagged as not converged.
build_row <- function(model_name, fit, n_pars) {
  if (is.null(fit))
    return(data.frame(clade = TARGET_CLADE, model = model_name, lambda = NA,
                      mu = NA, K = NA, loglik = NA, n_params = n_pars,
                      AIC = NA, converged = FALSE))
  data.frame(
    clade = TARGET_CLADE, model = model_name,
    lambda = if ("lambda" %in% names(fit)) fit$lambda else fit$lambda0,
    mu = if ("mu" %in% names(fit)) fit$mu else fit$mu0,
    K = if ("K" %in% names(fit)) fit$K else NA,
    loglik = fit$loglik, n_params = n_pars,
    AIC = -2 * fit$loglik + 2 * n_pars,
    converged = ifelse(is.null(fit$conv), NA, fit$conv == 0))
}

# Build the results table, with one row per model (the last argument is the number of estimated parameters).
results <- rbind(build_row("CR", cr_fit, 2),
                 build_row("DD_L", ddl_fit, 3),
                 build_row("DD_M", ddm_fit, 3))

# Remove the models that failed to fit, calculate the delta AIC and the Akaike weights, and sort the models from the best (lowest AIC) to the worst.
comparison <- results %>%
  filter(!is.na(AIC)) %>%
  mutate(delta_AIC = AIC - min(AIC),
         Akaike_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))) %>%
  arrange(AIC)

# Select the best-fitting model and check if its support over the second best one is strong (delta AIC > 2).
best <- comparison[1, ]
strong_support <- if (nrow(comparison) > 1) comparison$delta_AIC[2] > 2 else NA

# Create the interpretation message based on the best-fitting model.
msg <- switch(best$model,
              CR   = "No evidence of diversity-dependence - constant-rate model preferred.",
              DD_L = "Speciation rate appears to decline as the clade fills up.",
              DD_M = "Extinction rate appears to increase as the clade fills up.")

# If the support over the second best model is weak (delta AIC < 2), flag the result as inconclusive.
if (!is.na(strong_support) && !strong_support)
  msg <- paste(msg, "Support over the runner-up is weak (delta AIC < 2) - inconclusive.")

# If a diversity-dependent model is the best one, compare the estimated K with the number of known species:
# a K lower than three times the known diversity supports a genuine ceiling, while a far larger K means that the ceiling is poorly constrained.
if (best$model %in% c("DD_L", "DD_M") && !is.na(best$K))
  msg <- paste(msg, if (best$K / KNOWN_SPECIES < 3)
    "Estimated K is close to known diversity, supporting a genuine ceiling."
    else "Estimated K is far larger than known diversity - ceiling poorly constrained.")

# Add the interpretation message to the comparison table.
comparison$interpretation <- msg

# Save the comparison table of the selected clade as a CSV file.
write.csv(comparison,
          file.path(output_dir, paste0(TARGET_CLADE, "_dd_comparison.csv")),
          row.names = FALSE)

# Define the combined CSV file that holds the results of all the clades analysed so far.
combined_file <- file.path(output_dir, "ALL_clades_comparison.csv")

# If the combined file already exists, append the new results to the previous ones; otherwise, create it.
# Any earlier run of the selected clade is dropped, so that re-running a clade replaces its rows rather than duplicating them.
if (file.exists(combined_file)) {
  previous <- read.csv(combined_file, stringsAsFactors = FALSE)
  previous <- previous[previous$clade != TARGET_CLADE, ]
  combined <- rbind(previous, comparison)
} else {
  combined <- comparison
}

# Save the updated combined table.
write.csv(combined, combined_file, row.names = FALSE)

# Print the results of the selected clade as screen output.
cat("\n--- RESULTS:", TARGET_CLADE, "---\n")
print(comparison %>% select(model, lambda, mu, K, loglik, AIC, delta_AIC,
                            Akaike_weight, converged))
cat("\n", msg, "\n")
cat("\nSaved. Combined table now holds",
    length(unique(combined$clade)), "clade(s).\n")
