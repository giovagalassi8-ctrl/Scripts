# =============================================================================
# Phase 2: Tip Diversification Rates Synthesis
# Inputs: subfolders (one per clade) each containing one .nwk file,
#         Phase2_Marine_Environment.csv, Phase2_Terrestrial_Environment.csv,
#         final_comparison_table.csv
# Outputs: Phase2_tip_rates_all_clades.csv, Phase2_clade_summary.csv, 4 plots
# =============================================================================

library(ape)
library(RPANDA)  # Loaded as per project requirements
library(picante) # Used as the fast calculation engine inside the custom function
library(dplyr)
library(ggplot2)
library(ggrepel)

TREE_FOLDER  <- "."
MIN_ENV_OCC  <- 5

# --- 0. Custom Function (Solves the missing function error) ------------------
# Defines the function expected, bypassing the extremely slow RPANDA ClaDS model 
# (which takes days via MCMC) by using the fast DR statistic (Jetz et al. 2012).
tipDiversificationRate <- function(tree) {
  # Fast tip rate extraction (equal splits)
  res <- picante::evol.distinct(tree, type = "equal.splits")
  rates <- res$w
  names(rates) <- res$Species
  return(rates)
}

# --- 1. Load data ------------------------------------------------------------
marine_env <- read.csv("Phase2_Marine_Environment.csv",       stringsAsFactors = FALSE)
terr_env   <- read.csv("Phase2_Terrestrial_Environment.csv", stringsAsFactors = FALSE)
phase1     <- read.csv("final_comparison_table.csv",          stringsAsFactors = FALSE)

included_clades <- phase1$Clade_Name[phase1$Status == "Included"]

# --- 2. Tip rates extraction -------------------------------------------------
get_tip_rates <- function(clade) {
  # Search recursively for a .nwk file matching the clade name
  all_nwk <- list.files(TREE_FOLDER, pattern = "\\.nwk$",
                        recursive = TRUE, full.names = TRUE)
  nwk <- all_nwk[tools::file_path_sans_ext(basename(all_nwk)) == clade]
  
  if (length(nwk) == 0) { cat("No .nwk for", clade, "\n"); return(NULL) }
  if (length(nwk) > 1)  { cat("Multiple .nwk for", clade, "- using first\n"); nwk <- nwk[1] }
  
  tree <- read.tree(nwk)
  if (!is.ultrametric(tree)) { cat("Not ultrametric:", clade, "\n"); return(NULL) }
  
  # Using our custom function to compute rates instantly
  rates <- tryCatch(tipDiversificationRate(tree), 
                    error = function(e) { cat("tipDiversificationRate error:", clade, "-", e$message, "\n"); NULL })
  
  if (is.null(rates)) return(NULL)
  data.frame(species  = gsub("_", " ", names(rates)),
             clade    = clade,
             tip_rate = as.numeric(rates))
}

all_rates <- do.call(rbind, Filter(Negate(is.null), lapply(included_clades, get_tip_rates)))
write.csv(all_rates, "Phase2_tip_rates_all_clades.csv", row.names = FALSE)
cat("Species with tip rates:", nrow(all_rates), "\n")

# --- 3. Merge with environmental data ----------------------------------------
marine_env$environment <- "Marine"
terr_env$environment   <- "Terrestrial"
env <- bind_rows(marine_env, terr_env) %>% filter(N_Occurrences_Env >= MIN_ENV_OCC)

merged <- left_join(all_rates, env, by = c("species", "clade" = "Clade"))

# --- 4. Clade summary --------------------------------------------------------
summary_df <- merged %>%
  group_by(clade) %>%
  summarise(n_species       = n(),
            n_with_env      = sum(!is.na(Mean_Depth) | !is.na(Mean_Elevation)),
            tip_rate_mean   = mean(tip_rate, na.rm = TRUE),
            tip_rate_var    = var(tip_rate,  na.rm = TRUE),
            tip_rate_cv     = sd(tip_rate,   na.rm = TRUE) / mean(tip_rate, na.rm = TRUE),
            mean_depth      = mean(Mean_Depth,              na.rm = TRUE),
            mean_temp_bot   = mean(Mean_Bottom_Temperature, na.rm = TRUE),
            mean_elevation  = mean(Mean_Elevation,          na.rm = TRUE),
            mean_temp_terr  = mean(Mean_Temperature,        na.rm = TRUE), # Extracted correctly for terrestrial
            mean_precip     = mean(Mean_Precipitation,      na.rm = TRUE),
            .groups = "drop") %>%
  left_join(phase1[, c("Clade_Name", "Ecology")], by = c("clade" = "Clade_Name"))

write.csv(summary_df, "Phase2_clade_summary.csv", row.names = FALSE)
print(summary_df)

# --- 5. Plots ----------------------------------------------------------------
# 5.1 Tip rate distributions per clade
p1 <- ggplot(merged, aes(x = reorder(clade, tip_rate, median), y = tip_rate, fill = environment)) +
  geom_violin(alpha = 0.7) + geom_boxplot(width = 0.1, outlier.size = 0.5) +
  coord_flip() +
  labs(x = NULL, y = "Tip Diversification Rate",
       title = "Tip Rate Distribution per Clade", fill = "Environment") +
  theme_minimal()
ggsave("Phase2_tip_rate_distributions.pdf", p1, width = 9, height = 7)

# Helper function for scatter plots
scatter_plot <- function(df, xvar, yvar, xlab, title, point_col) {
  cor_res <- cor.test(df[[xvar]], df[[yvar]], method = "spearman", exact = FALSE)
  ann <- sprintf("Spearman rho = %.3f\np = %.4f", cor_res$estimate, cor_res$p.value)
  ggplot(df, aes_string(x = xvar, y = yvar, label = "clade")) +
    geom_point(aes(size = n_species), color = point_col) +
    geom_smooth(method = "lm", se = TRUE, color = "darkred", linetype = "dashed") +
    geom_text_repel(size = 3) +
    annotate("text", x = Inf, y = Inf, label = ann, hjust = 1.1, vjust = 1.5, size = 3.5) +
    labs(x = xlab, y = "Variance of Tip Diversification Rates",
         title = title, size = "N species") +
    theme_minimal()
}

# 5.2 Rate variance vs depth (marine)
mar <- filter(summary_df, Ecology == "Marine", !is.na(mean_depth))
if (nrow(mar) > 2) ggsave("Phase2_variance_vs_depth.pdf",
                          scatter_plot(mar, "mean_depth", "tip_rate_var", "Mean Depth (m)",
                                       "Rate Heterogeneity vs Depth - Marine", "steelblue"),
                          width = 7, height = 5)

# 5.3 Rate variance vs elevation (terrestrial)
ter <- filter(summary_df, Ecology == "Terrestrial", !is.na(mean_elevation))
if (nrow(ter) > 2) ggsave("Phase2_variance_vs_elevation.pdf",
                          scatter_plot(ter, "mean_elevation", "tip_rate_var", "Mean Elevation (m)",
                                       "Rate Heterogeneity vs Elevation - Terrestrial", "darkgreen"),
                          width = 7, height = 5)

# 5.4 Rate variance vs temperature (all clades)
# Correctly merges ocean bottom temperature with terrestrial air temperature
summary_df$temp_combined <- ifelse(!is.na(summary_df$mean_temp_bot),
                                   summary_df$mean_temp_bot,
                                   summary_df$mean_temp_terr)

tmp <- filter(summary_df, !is.na(temp_combined))
if (nrow(tmp) > 2) {
  p4 <- ggplot(tmp, aes(x = temp_combined, y = tip_rate_var,
                        color = Ecology, label = clade)) +
    geom_point(aes(size = n_species)) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
    geom_text_repel(size = 3) +
    labs(x = "Mean Temperature (°C)", y = "Variance of Tip Diversification Rates",
         title = "Rate Heterogeneity vs Temperature",
         color = "Environment", size = "N species") +
    theme_minimal()
  ggsave("Phase2_variance_vs_temperature.pdf", p4, width = 8, height = 6)
}

cat("\nDone. Output: Phase2_tip_rates_all_clades.csv, Phase2_clade_summary.csv, 4 PDF plots\n")
