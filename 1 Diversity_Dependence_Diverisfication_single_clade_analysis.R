# =============================================================================
# Diversity-Dependent Diversification (DDD) - single clade version
#
# WHAT THIS SCRIPT DOES
# Fits and compares three diversification models for ONE clade at a time:
#   CR   - Constant Rate: rates fixed over time (null model, 2 parameters)
#   DD_L - speciation DECLINES as lineages approach carrying capacity K
#   DD_M - extinction INCREASES as lineages approach carrying capacity K
# Models are compared with AIC and Akaike weights.
#
# HOW TO USE IT
#   1. Set the two parameters below: the clade name and its known species count.
#      The folder is found automatically by searching DATA_ROOT for a folder
#      matching the clade name, whatever phylum it sits under.
#   2. Run the script (use Run / Ctrl+Enter, not source(), so progress
#      messages appear as they happen).
#   3. When it finishes, change the two parameters and run again.
#
# Results accumulate: each run updates a combined CSV holding every clade
# analysed so far, so you never need to re-run earlier ones.
#
# OUTPUTS (in 04_DDD/)
#   <clade>_dd_fits.rds        - raw model fits for this clade
#   <clade>_dd_comparison.csv  - AIC comparison for this clade
#   ALL_clades_comparison.csv  - combined table, updated after every run

# Set the initial lamba value for the analysis (change accordingly).
#initial_lambda <- 0.3
# Set the initial mu value for the analysis (change accordingly).
#initial_mu     <- 0.1
# =============================================================================


library(ape)
library(DDD)
library(dplyr)

# Set the clade to analyse (change accordingly).
# Expected structure of the repository: there are several numbererd sub-folders of which one, usually the fist, contains the data (in this case, it is called 00_Data).
# The data folder contains other sub-folders, one per clade, containing all the data for each group, including a treefile in nwk format for that group.
# The name selected below must match the clade data folder name exactly.
TARGET_CLADE <- "CLADE_NAME" 

# Set the total number of total accepted extant species (taken from online databases, like WoRMS, or from literature).
KNOWN_SPECIES <- 1000


# Set the data folder (change accordingly).
DATA_ROOT <- "DATA_FOLDER_NAME"
# Set the output directory name (change accordingly).
output_dir     <- "OUTPUT_DIRECTORY_NAME"

# Create the output directory if it does not exixst.
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Searches every subfolder in the DATA_ROOT folder and keeps the one whose name matches the clade.
all_dirs <- list.dirs(DATA_ROOT, recursive = TRUE, full.names = TRUE)
clade_folder <- all_dirs[basename(all_dirs) == TARGET_CLADE]

# Stop the process if there is not a match between the name selected and the name of any of the subfolders.
if (length(clade_folder) == 0)
  stop("No folder named '", TARGET_CLADE, "' found under ", DATA_ROOT)
# Stop the process if there are multiple subfolders with the same matching name.
if (length(clade_folder) > 1)
  stop("Multiple folders named '", TARGET_CLADE, "' found - resolve the ambiguity.")


# Look for the newick format treefile into the matching-name folder for the selected clade.
# Stop if the newick treefile has not found.
nwk_files <- list.files(clade_folder, pattern = "\\.nwk$", full.names = TRUE)
if (length(nwk_files) == 0)
  stop("No .nwk file found in ", clade_folder)
# Read the treefile for further analyses and check if it contains more than one tree (in this case it uses only the first one) and if the tree is ultrametric.
tree <- read.tree(nwk_files[1])
if (inherits(tree, "multiPhylo")) tree <- tree[[1]]
if (!is.ultrametric(tree))
  stop("Tree is not ultrametric - DDD requires an ultrametric tree.")

# Create an object containing the number of tips of the tree.
n_tips <- Ntip(tree)
# Calculate the number of missing species comparing the number of species in the tree with the number of accepted extant species settled before.
missnumspec <- max(0, KNOWN_SPECIES - n_tips)
# Create an object containing the branching times of the tree.
brts <- branching.times(tree)
# Returns a summary of some tree's statistics as screen output. 
# In particular, the ratio is returned: given the large ratio value, the script may take much longer, as it is more complex to manage the model (especially for values greater than 12).
cat("Tips:", n_tips, "| Known species:", KNOWN_SPECIES, "| Missing:", missnumspec, "| Ratio:", round(missnumspec / n_tips, 1), "\n")

# Check if the tree has more tips than the number of known species, and stop the proccess if necessary.
if (n_tips > KNOWN_SPECIES)
  warning("Tree has MORE tips than the stated known species count.")
# Check if the ratio is higher than 12 and in case it warns you that the process could require a lot of time.
if (missnumspec / n_tips > 12)
  cat("WARNING: missing species greatly outnumber sampled tips.",
      "This fit may be slow and K poorly constrained.\n")

# --- Derive starting values from the tree itself ------------------------------
# A fixed starting guess (e.g. lambda=0.3 for every clade) can fail on large
# or fast-diversifying trees: the likelihood at that starting point can
# underflow to -Inf, or the optimiser can wander into a degenerate loop
# (as happened with Dugesiidae, 215 tips). Deriving starting values from the
# tree itself keeps the optimiser in a plausible region for THIS clade's
# actual scale, without changing what the true maximum likelihood is - it
# only changes where the search begins.
# Set the initial lamba value for the analysis (change accordingly).
#initial_lambda <- 0.3
# Set the initial mu value for the analysis (change accordingly).
#initial_mu     <- 0.1

tree_age       <- max(brts)
initial_lambda <- log(n_tips) / tree_age
# ^ rough net diversification rate estimate: ln(N tips) / crown age
initial_mu     <- initial_lambda * 0.3
# ^ extinction assumed at a fraction of speciation as a starting guess only

cat("Derived starting values: lambda =", round(initial_lambda, 4),
    "| mu =", round(initial_mu, 4), "\n")


# Fit the three models
#
cond      <- 1                          # condition on clade survival
initial_K <- (n_tips + missnumspec) * 2  # start K above total diversity

cat("\nFitting CR...\n")
cr_fit <- tryCatch(
  bd_ML(brts = brts, initparsopt = c(initial_lambda, initial_mu),
        idparsopt = c(1, 2), missnumspec = missnumspec, cond = cond),
  error = function(e) { cat("  failed:", e$message, "\n"); NULL })

cat("Fitting DD_L (diversity-dependent speciation)...\n")
ddl_fit <- tryCatch(
  dd_ML(brts = brts, initparsopt = c(initial_lambda, initial_mu, initial_K),
        idparsopt = c(1, 2, 3), missnumspec = missnumspec, cond = cond,
        ddmodel = 1),
  error = function(e) { cat("  failed:", e$message, "\n"); NULL })

cat("Fitting DD_M (diversity-dependent extinction)...\n")
ddm_fit <- tryCatch(
  dd_ML(brts = brts, initparsopt = c(initial_lambda, initial_mu, initial_K),
        idparsopt = c(1, 2, 3), missnumspec = missnumspec, cond = cond,
        ddmodel = 3),
  error = function(e) { cat("  failed:", e$message, "\n"); NULL })

saveRDS(list(clade = TARGET_CLADE, n_tips = n_tips, missnumspec = missnumspec,
             known_species = KNOWN_SPECIES,
             CR = cr_fit, DD_L = ddl_fit, DD_M = ddm_fit),
        file.path(output_dir, paste0(TARGET_CLADE, "_dd_fits.rds")))

# --- Build the comparison table ------------------------------------------------

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

results <- rbind(build_row("CR", cr_fit, 2),
                 build_row("DD_L", ddl_fit, 3),
                 build_row("DD_M", ddm_fit, 3))

comparison <- results %>%
  filter(!is.na(AIC)) %>%
  mutate(delta_AIC = AIC - min(AIC),
         Akaike_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC))) %>%
  arrange(AIC)

# --- Interpret ------------------------------------------------------------------

best <- comparison[1, ]
strong_support <- if (nrow(comparison) > 1) comparison$delta_AIC[2] > 2 else NA

msg <- switch(best$model,
              CR   = "No evidence of diversity-dependence - constant-rate model preferred.",
              DD_L = "Speciation rate appears to decline as the clade fills up.",
              DD_M = "Extinction rate appears to increase as the clade fills up.")

if (!is.na(strong_support) && !strong_support)
  msg <- paste(msg, "Support over the runner-up is weak (delta AIC < 2) - inconclusive.")

if (best$model %in% c("DD_L", "DD_M") && !is.na(best$K))
  msg <- paste(msg, if (best$K / KNOWN_SPECIES < 3)
    "Estimated K is close to known diversity, supporting a genuine ceiling."
    else "Estimated K is far larger than known diversity - ceiling poorly constrained.")

comparison$interpretation <- msg

write.csv(comparison,
          file.path(output_dir, paste0(TARGET_CLADE, "_dd_comparison.csv")),
          row.names = FALSE)

# --- Append to the combined table across all clades run so far ------------------

combined_file <- file.path(output_dir, "ALL_clades_comparison.csv")

if (file.exists(combined_file)) {
  previous <- read.csv(combined_file, stringsAsFactors = FALSE)
  previous <- previous[previous$clade != TARGET_CLADE, ]
  # ^ drops any earlier run of this clade, so re-running replaces rather
  #   than duplicating its rows
  combined <- rbind(previous, comparison)
} else {
  combined <- comparison
}

write.csv(combined, combined_file, row.names = FALSE)

# --- Report ---------------------------------------------------------------------

cat("\n--- RESULTS:", TARGET_CLADE, "---\n")
print(comparison %>% select(model, lambda, mu, K, loglik, AIC, delta_AIC,
                            Akaike_weight, converged))
cat("\n", msg, "\n")
cat("\nSaved. Combined table now holds",
    length(unique(combined$clade)), "clade(s).\n")
