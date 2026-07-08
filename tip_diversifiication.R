#!/usr/bin/env Rscript

# This script synthesizes tip diversification rates for various animal clades.
# It bypasses computationally heavy models by calculating the fast DR statistic using the 'picante' package.
# It then merges these evolutionary rates with the previously extracted marine and terrestrial environmental data, computes clade-level summaries, and 
# generates comparative plots (violins and scatter plots).
# It requires the csv file containing the environmental parameters obtained by 'terrestrial_environment_extraction.R' and 'marine_environment_extraction.R' scripts.

library(ape)
library(RPANDA)
library(picante)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Define the root directory path where the script will search for the phylogenetic tree files.
TREE_FOLDER  <- "."
# Define the minimum threshold of environmental occurrences required to retain a species in the final analysis.
MIN_ENV_OCC  <- 5

# Define a custom function ('tipDiversificationRate') that takes a phylogenetic 'tree' object as input.
tipDiversificationRate <- function(tree) {
  # Calculate the evolutionary distinctiveness (DR statistic) for each tip using the "equal.splits" method from 'picante'.
  res <- picante::evol.distinct(tree, type = "equal.splits")
  # Extract the numeric vector of calculated rates (the 'w' column from the result).
  rates <- res$w
  # Assign the species names (from the 'Species' column) as the names attribute of the 'rates' vector.
  names(rates) <- res$Species
  # Return the named numeric vector containing the tip diversification rates.
  return(rates)
}

# Read the marine environmental data CSV file into a dataframe.
marine_env <- read.csv("Marine_Environment.csv", stringsAsFactors = FALSE)
# Read the terrestrial environmental data CSV file into a dataframe.
terr_env <- read.csv("Terrestrial_Environment.csv", stringsAsFactors = FALSE)
# Read the metadata CSV file obtained by ... script (which contains the inclusion status of clades) into a dataframe
comparison <- read.csv("final_comparison_table.csv", stringsAsFactors = FALSE)

# Subset the 'Clade_Name' column keeping only the rows where the 'Status' column equals "Included".
included_clades <- comparison$Clade_Name[comparison$Status == "Included"]

# Define a function (get_tip_rates).
get_tip_rates <- function(clade) {
  # Search recursively for a .nwk file matching the clade name.
  all_nwk <- list.files(TREE_FOLDER, pattern = "\\.nwk$",
                        recursive = TRUE, full.names = TRUE)
  # Filter 'all_nwk' to find the specific file where the base name (without extension) exactly matches the 'clade' name.
  nwk <- all_nwk[tools::file_path_sans_ext(basename(all_nwk)) == clade]

  # Check if the length of 'nwk' is 0 (no file found); if so, print a warning and return NULL to skip this clade.
  if (length(nwk) == 0) { cat("No .nwk for", clade, "\n"); return(NULL) }
  # Check if more than 1 file was found; if so, print a warning and keep only the first file in the vector. 
  if (length(nwk) > 1)  { cat("Multiple .nwk for", clade, "- using first\n"); nwk <- nwk[1] }
  
  tree <- read.tree(nwk)
  # Check if the tree is NOT ultrametric (tips don't align at time 0).
  if (!is.ultrametric(tree)) { cat("Not ultrametric:", clade, "\n"); return(NULL) }
  
  # Using our custom function to compute rates instantly.
  rates <- tryCatch(tipDiversificationRate(tree), 
                    error = function(e) { cat("tipDiversificationRate error:", clade, "-", e$message, "\n"); NULL })
  # Check if 'rates' is NULL (meaning the calculation failed).
  if (is.null(rates)) return(NULL)
  # Create and return a new dataframe containing the extracted species, clade names, and their rates.
  data.frame(
    # All underscores ('_') are replaced with spaces (' ') in the names of the rates vector.
    species = gsub("_", " ", names(rates)),
    # Create the 'clade' column.
    clade = clade,
    # Create the 'tip_rate' column by converting the rates vector into a pure numeric format.
    tip_rate = as.numeric(rates))
}

# Apply 'get_tip_rates' to all 'included_clades', remove NULL results using 'Filter' and 'Negate', and row-bind the rest into 'all_rates'.
all_rates <- do.call(rbind, Filter(Negate(is.null), lapply(included_clades, get_tip_rates)))
# Save the tip rates file.
write.csv(all_rates, "tip_rates_all_clades.csv", row.names = FALSE)
cat("Species with tip rates:", nrow(all_rates), "\n")

# Merge with environmental data.
# Add a new column named 'environment' to the marine and terrestrial dataframe and fill it with the setted string.
marine_env$environment <- "Marine"
terr_env$environment   <- "Terrestrial"
# Row-bind the marine and terrestrial dataframes together, and keep only the rows where the number of environmental occurrences is greater than or equal to MIN_ENV_OCC.
env <- bind_rows(marine_env, terr_env) %>% filter(N_Occurrences_Env >= MIN_ENV_OCC)

# Perform a left join to attach the filtered environmental data ('env') to the evolutionary rates ('all_rates') matching by species and clade.
merged <- left_join(all_rates, env, by = c("species", "clade" = "Clade"))

# Create the summary file by starting with the 'merged' dataframe.
summary_df <- merged %>%
  # Group the data by the 'clade' column.
  group_by(clade) %>%
  summarise(
            # Count the total number of species in the current clade group.
            n_species = n(),
            # Count how many species in the clade have valid environmental data.
            n_with_env = sum(!is.na(Mean_Depth) | !is.na(Mean_Elevation)),
            # Calculate the mean of the tip diversification rates, ignoring NA values.
            tip_rate_mean   = mean(tip_rate, na.rm = TRUE),
            # Calculate the variance of the tip diversification rates, ignoring NA values.
            tip_rate_var    = var(tip_rate,  na.rm = TRUE),
            # Calculate the Coefficient of Variation (standard deviation divided by mean) for the tip rates.
            tip_rate_cv     = sd(tip_rate,   na.rm = TRUE) / mean(tip_rate, na.rm = TRUE),
            # Calculate the mean of the depth values for the clade, ignoring NA values.
            mean_depth      = mean(Mean_Depth,              na.rm = TRUE),
            # Calculate the mean of the bottom temperature values for the clade, ignoring NA values.
            mean_temp_bot   = mean(Mean_Bottom_Temperature, na.rm = TRUE),
            # Calculate the mean of the elevation values for the clade, ignoring NA values.
            mean_elevation  = mean(Mean_Elevation,          na.rm = TRUE),
            # Calculate the mean of the terrestrial temperature values for the clade, ignoring NA values.
            mean_temp_terr  = mean(Mean_Temperature,        na.rm = TRUE),
            # Calculate the mean of the precipitation values for the clade, ignoring NA values.
            mean_precip     = mean(Mean_Precipitation,      na.rm = TRUE),
            # Drop the grouping structure of the resulting dataframe to prevent issues in downstream joins 
            .groups = "drop") %>%
  left_join(comparison[, c("Clade_Name", "Ecology")], by = c("clade" = "Clade_Name"))

# Save the final statistical summary dataframe.
write.csv(summary_df, "clade_summary.csv", row.names = FALSE)
print(summary_df)

# Tip rate distributions per clade (Violin Plot).
# Initialize a ggplot object using the 'merged' dataframe, mapping clade (reordered by median tip_rate) to x, tip_rate to y, and environment to fill color.
p1 <- ggplot(merged, aes(x = reorder(clade, tip_rate, median), y = tip_rate, fill = environment)) +
  # Add violin plots.
  geom_violin(alpha = 0.7) + 
  # Add narrow boxplots over the violins to show quartiles and outliers.
  geom_boxplot(width = 0.1, outlier.size = 0.5) +
  # Swap the x and y axes to make the long clade names readable horizontally.
  coord_flip() +
  labs(x = NULL,
       y = "Tip Diversification Rate",
       title = "Tip Rate Distribution per Clade",
       fill = "Environment") +
  # Apply a clean, minimalist theme to the plot.
  theme_minimal()
# Save the plot.
ggsave("tip_rate_distributions.pdf", p1, width = 9, height = 7)

# Define a helper function (scatter_plot) to standardise the creation of correlation scatter plots
scatter_plot <- function(df, xvar, yvar, xlab, title, point_col) {
  # Perform a Spearman rank correlation test between the x and y variables, disabling exact p-value calculation to avoid warnings with ties.
  cor_res <- cor.test(df[[xvar]], df[[yvar]], method = "spearman", exact = FALSE)
  # Format a text string containing the Spearman rho coefficient and the p-value for the plot annotation.
  ann <- sprintf("Spearman rho = %.3f\np = %.4f", cor_res$estimate, cor_res$p.value)

  # Creates the plot.
  ggplot(df, aes_string(x = xvar, 
                        y = yvar, 
                        label = "clade")) +
    # Add points to the scatter plot.
    geom_point(aes(size = n_species),
               color = point_col) +
    # Add a linear regression trendline without the standard error shading.
    geom_smooth(method = "lm", 
                se = FALSE,
                color = "darkred",
                linetype = "dashed") +
    # Add non-overlapping text labels to the points using ggrepel.
    geom_text_repel(size = 3) +
    # Add a custom text annotation.
    annotate("text",
             x = Inf,
             y = Inf, 
             label = ann,
             hjust = 1.1, 
             vjust = 1.5,
             size = 3.5) +
    labs(x = xlab,
         y = "Variance of Tip Diversification Rates",
         title = title,
         size = "N species") +
    # Apply a clean, minimalist theme to the plot.
    theme_minimal()
}

# Rate variance vs. Mean Depth (Marine only).
# Filter the summary dataframe to keep only marine clades that have valid mean depth values.
mar <- filter(summary_df, Ecology == "Marine", !is.na(mean_depth))
# Check if there are more than 2 marine clades (correlation requires at least 3 points).
if (nrow(mar) > 2) ggsave("variance_vs_depth.pdf",
                          scatter_plot(mar, "mean_depth", "tip_rate_var", "Mean Depth (m)",
                                       "Rate Heterogeneity vs Depth - Marine", "steelblue"),
                          width = 7, height = 5)

# Rate variance vs. Mean Elevation (Terrestrial only).
# Filter the summary dataframe to keep only terrestrial clades that have valid mean elevation values.
ter <- filter(summary_df, Ecology == "Terrestrial", !is.na(mean_elevation))
# Check if there are more than 2 terrestrial clades (correlation requires at least 3 points).
if (nrow(ter) > 2) ggsave("variance_vs_elevation.pdf",
                          scatter_plot(ter, "mean_elevation", "tip_rate_var", "Mean Elevation (m)",
                                       "Rate Heterogeneity vs Elevation - Terrestrial", "darkgreen"),
                          width = 7, height = 5)

# Rate variance vs. Temperature (All clades).
# Correctly merges ocean bottom temperature with terrestrial air temperature.
summary_df$temp_combined <- ifelse(!is.na(summary_df$mean_temp_bot),
                                   summary_df$mean_temp_bot,
                                   summary_df$mean_temp_terr)
# Filter the summary dataframe to keep only rows with a valid combined temperature value.
tmp <- filter(summary_df, !is.na(temp_combined))
# Check if there are more than 2 valid clades left for the temperature analysis.
if (nrow(tmp) > 2) {
  # Creates the plot.
  p4 <- ggplot(tmp, aes(x = temp_combined,
                        y = tip_rate_var,
                        color = Ecology,
                        label = clade)) +
    # Add points, mapping point size to the number of species.
    geom_point(aes(size = n_species)) +
    # Add linear trendlines for each ecology group without standard error shading.
    geom_smooth(method = "lm",
                se = FALSE,
                linetype = "dashed") +
    # Add non-overlapping text labels to the points.
    geom_text_repel(size = 3) +
    labs(x = "Mean Temperature (°C)",
         y = "Variance of Tip Diversification Rates",
         title = "Rate Heterogeneity vs Temperature",
         color = "Environment",
         size = "N species") +
    # Apply a clean, minimalist theme to the plot.
    theme_minimal()
  
  # Save the temperature plot.
  ggsave("variance_vs_temperature.pdf", p4, width = 7, height = 5)
}
