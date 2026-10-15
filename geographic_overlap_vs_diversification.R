#!/usr/bin/env Rscript

# This script tests whether the geographic overlap between the species of a clade is related to its diversification dynamics,
# comparing the mean geographic overlap of each clade with its Gamma and DeltaR statistics.
# It produces:
#   1. a diagnostic plot of all the clades, highlighting the ones excluded from the analysis;
#   2. the main plots (Gamma and DeltaR) testing the general hypothesis on all the included clades, with the results of the Spearman correlation;
#   3. supplementary plots faceted by ecology (marine, terrestrial, freshwater), which decompose the trends to reveal confounding spatial scales
#      (justifying a further analysis on depth).
# Clades are excluded from the analysis if they are parasites, if their tree has less than 15 tips or if their sampling fraction is lower than 0.05.
# It requires a CSV file containing one row per clade, with the clade name, the mean geographic overlap (obtained with the 'geographic_overlap.R' script),
# the number of tips, the sampling fraction and the Gamma and DeltaR statistics (obtained with the 'ltt_plot_and_gamma.R' script).
# The table of the included clades is saved as a CSV file.


library(dplyr)
library(ggplot2)
library(ggrepel)

# Check if the input file exists, otherwise stop the process (change with the correct file path, here and in the following command).
if (!file.exists("analyses/Geographic_overlap/All_data.csv")) {
  stop("Input file 'All_data.csv' not found.")
}

# Read the CSV file containing the data of all the clades.
clade_data <- read.csv("analyses/Geographic_overlap/All_data.csv", stringsAsFactors = FALSE)
# Convert the column names into syntactically valid R names (e.g. spaces are replaced by dots).
colnames(clade_data) <- make.names(colnames(clade_data))

# Assign each clade to an ecological category and define its inclusion status in the analysis (change the clade names accordingly).
# Clades are excluded if they are parasites, if their tree has less than 15 tips (low statistical power) or if their sampling fraction is lower than 0.05.
clade_data <- clade_data %>%
  mutate(
    Ecology = case_when(
      Clade_Name %in% c("Carditidae", "Conidae", "Siboglinidae", "Strombidae",
                        "Terebratulina", "Terebridae", "Veneridae", "Monocelididae", "Lingula") ~ "Marine",
      Clade_Name %in% c("Arion", "Clausiliidae", "Haemadipsidae", "Hormogastridae",
                        "Lumbricidae", "Theba") ~ "Terrestrial",
      Clade_Name %in% c("Dugesiidae", "Hirudo", "Unionidae") ~ "Freshwater",
      Clade_Name %in% c("Polystomatidae", "Taeniidae") ~ "Parasite",
      TRUE ~ "Other"
    ),
    Status = case_when(
      Ecology == "Parasite" ~ "Excluded (Parasite)",
      n.tips < 15 ~ "Excluded (Low Power)",
      sampling.fraction < 0.05 ~ "Excluded (Low Sampling)",
      TRUE ~ "Included"
    )
  )

# Create the diagnostic plot of all the clades, in which crosses indicate the excluded ones.
plot_diagnostic <- ggplot(clade_data, aes(x = Mean_Geographic_Overlap, y = gamma)) +
  # Draw the points, colored by ecology, sized by the number of tips and shaped by inclusion status.
  geom_point(aes(color = Ecology, size = n.tips, shape = Status == "Included"), alpha = 0.8) +
  # Add the clade names as labels, avoiding overlaps.
  geom_text_repel(aes(label = Clade_Name, color = Ecology),
                  size = 3, fontface = "bold", max.overlaps = 30, show.legend = FALSE) +
  # Use circles for the included clades and crosses for the excluded ones.
  scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 4),
                     labels = c("TRUE" = "Valid for Analysis", "FALSE" = "Excluded")) +
  # Assign a specific color to each ecological category (change accordingly).
  scale_color_manual(values = c("Marine"="#1f78b4", "Terrestrial"="#33a02c",
                                "Freshwater"="#ff7f00", "Parasite"="grey50")) +
  # Set the title, subtitle, axis and legend labels.
  labs(
    x = "Mean Geographic Overlap (Jaccard Index)", y = "Gamma Statistic",
    title = "Diagnostic Overview",
    subtitle = "(Crosses indicate clades excluded due to parasite lifestyle or low statistical power)",
    shape = "Data Quality",
    color = "Environment",
    size = "Number of Tips"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# Show the plot.
plot(plot_diagnostic)

# Keep only the clades included in the analysis.
clade_data_filtered <- clade_data %>% filter(Status == "Included")

# Calculate the Spearman correlation between the mean geographic overlap and the Gamma statistic, and between the mean geographic overlap and DeltaR.
cor_gamma <- cor.test(clade_data_filtered$Mean_Geographic_Overlap, clade_data_filtered$gamma, method = "spearman", exact = FALSE)
cor_deltar <- cor.test(clade_data_filtered$Mean_Geographic_Overlap, clade_data_filtered$deltaR, method = "spearman", exact = FALSE)

# Create the text labels with the correlation coefficient (rho) and the p-value, to be shown in the main plots.
stats_gamma_text <- paste0("Spearman rho = ", round(cor_gamma$estimate, 3), "\np-value = ", round(cor_gamma$p.value, 4))
stats_deltar_text <- paste0("Spearman rho = ", round(cor_deltar$estimate, 3), "\np-value = ", round(cor_deltar$p.value, 4))

# Create the main plot of the mean geographic overlap vs the Gamma statistic.
plot_main_gamma <- ggplot(clade_data_filtered, aes(x = Mean_Geographic_Overlap, y = gamma)) +
  # Add a linear regression line with its confidence interval.
  geom_smooth(method = "lm", se = TRUE, color = "darkgrey", linetype = "dashed", alpha = 0.1) +
  # Draw the points, colored by ecology and sized by the number of tips.
  geom_point(aes(color = Ecology, size = n.tips), alpha = 0.9) +
  # Add the clade names as labels, avoiding overlaps.
  geom_text_repel(aes(label = Clade_Name, color = Ecology), size = 3.5, fontface = "bold", show.legend = FALSE) +
  # Assign a specific color to each ecological category (change accordingly).
  scale_color_manual(values = c("Marine"="#1f78b4", "Terrestrial"="#33a02c", "Freshwater"="#ff7f00")) +
  # Add the results of the Spearman correlation in the upper right part of the plot.
  annotate("text", x = max(clade_data_filtered$Mean_Geographic_Overlap)*0.75, y = max(clade_data_filtered$gamma)*0.9,
           label = stats_gamma_text, size = 4, fontface = "italic", hjust = 0) +
  # Set the title, axis and legend labels.
  labs(title = "Geographic Overlap vs Gamma",
       x = "Mean Geographic Overlap",
       y = "Gamma Statistic",
       size = "Number of Tips",
       color = "Environment") +
  theme_minimal()

# Show the plot.
plot(plot_main_gamma)

# Create the main plot of the mean geographic overlap vs DeltaR (the layers are the same as in the previous plot).
plot_main_deltar <- ggplot(clade_data_filtered, aes(x = Mean_Geographic_Overlap, y = deltaR)) +
  geom_smooth(method = "lm", se = TRUE, color = "darkgrey", linetype = "dashed", alpha = 0.1) +
  geom_point(aes(color = Ecology, size = n.tips), alpha = 0.9) +
  geom_text_repel(aes(label = Clade_Name, color = Ecology), size = 3.5, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c("Marine"="#1f78b4", "Terrestrial"="#33a02c", "Freshwater"="#ff7f00")) +
  annotate("text", x = max(clade_data_filtered$Mean_Geographic_Overlap)*0.75, y = max(clade_data_filtered$deltaR)*0.9,
           label = stats_deltar_text, size = 4, fontface = "italic", hjust = 0) +
  labs(title = "Geographic Overlap vs DeltaR",
       x = "Mean Geographic Overlap",
       y = "DeltaR",
       size = "Number of Tips",
       color = "Environment") +
  theme_minimal()

# Show the plot.
plot(plot_main_deltar)

# The following supplementary plots are faceted by ecology.
# P-values are omitted here because the sample size per facet (N = 3-5) is too small for a robust correlation:
# the regression lines are shown only for a qualitative visualization of the trends.

# Create the supplementary plot of the mean geographic overlap vs the Gamma statistic, faceted by ecology.
plot_supp_gamma <- ggplot(clade_data_filtered, aes(x = Mean_Geographic_Overlap, y = gamma)) +
  # Add a linear regression line without confidence interval.
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_point(aes(color = Ecology, size = n.tips), alpha = 0.9) +
  geom_text_repel(aes(label = Clade_Name, color = Ecology), size = 3.5, fontface = "bold", show.legend = FALSE) +
  # Split the plot into one panel per ecological category (fixed scales highlight the different spatial domains).
  facet_wrap(~ Ecology, scales = "fixed") +
  scale_color_manual(values = c("Marine"="#1f78b4", "Terrestrial"="#33a02c", "Freshwater"="#ff7f00")) +
  labs(title = "Overlap vs Gamma by Environment",
       x = "Mean Geographic Overlap",
       y = "Gamma Statistic",
       size = "Number of Tips") +
  theme_bw(base_size = 14) +
  # Hide the legend, since the facets are self-explanatory.
  theme(legend.position = "none")

# Show the plot.
plot(plot_supp_gamma)

# Create the supplementary plot of the mean geographic overlap vs DeltaR, faceted by ecology (the layers are the same as in the previous plot).
plot_supp_deltar <- ggplot(clade_data_filtered, aes(x = Mean_Geographic_Overlap, y = deltaR)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_point(aes(color = Ecology, size = n.tips), alpha = 0.9) +
  geom_text_repel(aes(label = Clade_Name, color = Ecology), size = 3.5, fontface = "bold", show.legend = FALSE) +
  facet_wrap(~ Ecology, scales = "fixed") +
  scale_color_manual(values = c("Marine"="#1f78b4", "Terrestrial"="#33a02c", "Freshwater"="#ff7f00")) +
  labs(title = "Overlap vs DeltaR by Environment",
       x = "Mean Geographic Overlap", y = "DeltaR", size = "Number of Tips") +
  theme_bw(base_size = 14) + theme(legend.position = "none")

# Show the plot.
plot(plot_supp_deltar)

# Save the table of the included clades as a CSV file (change the file name as desired).
write.csv(clade_data_filtered, "final_comparison_table.csv", row.names = FALSE)
