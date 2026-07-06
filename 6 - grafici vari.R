# =============================================================================
# Script: 04_compare_clades_final.R
# Purpose: Comprehensive macroevolutionary analysis pipeline.
# Outputs: 
#   1. Diagnostic plot (all 20 raw clades).
#   2. Main aggregate plots (Gamma & DeltaR) testing the universal hypothesis.
#   3. Supplementary faceted plots decomposing trends by Ecology to reveal
#      confounding spatial scales and justify Phase 2 (depth analysis).
# =============================================================================

library(dplyr)
library(ggplot2)
library(ggrepel)

# --- 1. LOAD AND PREPROCESS DATA ---------------------------------------------
if (!file.exists("All_data.csv")) {
  stop("Input file 'All_data.csv' not found.")
}

clade_data <- read.csv("All_data.csv", stringsAsFactors = FALSE)
colnames(clade_data) <- make.names(colnames(clade_data))

# --- 2. ECOLOGICAL CLASSIFICATION & FILTER STATUS ----------------------------
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

# --- 3. PLOT 0: DIAGNOSTIC RAW DATA ------------------------------------------
plot_diagnostic <- ggplot(clade_data, aes(x = Mean_Geographic_Overlap, y = gamma)) +
  geom_point(aes(color = Ecology, size = n.tips, shape = Status == "Included"), alpha = 0.8) +
  geom_text_repel(aes(label = Clade_Name, color = Ecology), 
                  size = 3, fontface = "bold", max.overlaps = 30, show.legend = FALSE) +
  scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 4), 
                     labels = c("TRUE" = "Valid for Analysis", "FALSE" = "Excluded")) +
  scale_color_manual(values = c("Marine"="#1f78b4", "Terrestrial"="#33a02c", 
                                "Freshwater"="#ff7f00", "Parasite"="grey50")) +
  labs(
    x = "Mean Geographic Overlap (Jaccard Index)", y = "Gamma Statistic (Raw)",
    title = "Plot 0: Diagnostic Overview (All 20 Clades)",
    subtitle = "Crosses indicate clades excluded due to parasite lifestyle or low statistical power",
    shape = "Data Quality", color = "Environment", size = "Number of Tips"
  ) +
  theme_minimal() + theme(legend.position = "right")

plot(plot_diagnostic)
ggsave("PLOT_0_Diagnostic.png", plot = plot_diagnostic, width = 10, height = 7, dpi = 300, bg = "white")

# --- 4. APPLY FILTERING & GLOBAL STATISTICS ----------------------------------
clade_data_filtered <- clade_data %>% filter(Status == "Included")

cor_gamma <- cor.test(clade_data_filtered$Mean_Geographic_Overlap, clade_data_filtered$gamma, method = "spearman", exact = FALSE)
cor_deltar <- cor.test(clade_data_filtered$Mean_Geographic_Overlap, clade_data_filtered$deltaR, method = "spearman", exact = FALSE)

stats_gamma_text <- paste0("Spearman rho = ", round(cor_gamma$estimate, 3), "\np-value = ", round(cor_gamma$p.value, 4))
stats_deltar_text <- paste0("Spearman rho = ", round(cor_deltar$estimate, 3), "\np-value = ", round(cor_deltar$p.value, 4))

# --- 5. MAIN PLOTS: AGGREGATED HYPOTHESIS TEST -------------------------------
# 5A: Gamma
plot_main_gamma <- ggplot(clade_data_filtered, aes(x = Mean_Geographic_Overlap, y = gamma)) +
  geom_smooth(method = "lm", se = TRUE, color = "darkgrey", linetype = "dashed", alpha = 0.1) +
  geom_point(aes(color = Ecology, size = n.tips), alpha = 0.9) +
  geom_text_repel(aes(label = Clade_Name, color = Ecology), size = 3.5, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c("Marine"="#1f78b4", "Terrestrial"="#33a02c", "Freshwater"="#ff7f00")) +
  annotate("text", x = max(clade_data_filtered$Mean_Geographic_Overlap)*0.75, y = max(clade_data_filtered$gamma)*0.9, 
           label = stats_gamma_text, size = 4, fontface = "italic", hjust = 0) +
  labs(title = "Main Fig 1: Geographic Overlap vs Gamma", subtitle = "Testing the universal hypothesis (All valid clades aggregated)",
       x = "Mean Geographic Overlap", y = "Gamma Statistic", size = "Number of Tips", color = "Environment") +
  theme_minimal()

plot(plot_main_gamma)
ggsave("PLOT_1_Main_Gamma.png", plot = plot_main_gamma, width = 8.5, height = 6, dpi = 300, bg = "white")

# 5B: DeltaR
plot_main_deltar <- ggplot(clade_data_filtered, aes(x = Mean_Geographic_Overlap, y = deltaR)) +
  geom_smooth(method = "lm", se = TRUE, color = "darkgrey", linetype = "dashed", alpha = 0.1) +
  geom_point(aes(color = Ecology, size = n.tips), alpha = 0.9) +
  geom_text_repel(aes(label = Clade_Name, color = Ecology), size = 3.5, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c("Marine"="#1f78b4", "Terrestrial"="#33a02c", "Freshwater"="#ff7f00")) +
  annotate("text", x = max(clade_data_filtered$Mean_Geographic_Overlap)*0.75, y = max(clade_data_filtered$deltaR)*0.9, 
           label = stats_deltar_text, size = 4, fontface = "italic", hjust = 0) +
  labs(title = "Main Fig 2: Geographic Overlap vs DeltaR", subtitle = "Testing the universal hypothesis (All valid clades aggregated)",
       x = "Mean Geographic Overlap", y = "DeltaR", size = "Number of Tips", color = "Environment") +
  theme_minimal()

plot(plot_main_deltar)
ggsave("PLOT_2_Main_DeltaR.png", plot = plot_main_deltar, width = 8.5, height = 6, dpi = 300, bg = "white")

# --- 6. SUPPLEMENTARY PLOTS: FACETED BY ECOLOGY ------------------------------
# We omit p-values here because the sample size per facet (N=3-5) is too small 
# for robust correlation. We show the regression lines only for qualitative trend visualization.

# 6A: Faceted Gamma
plot_supp_gamma <- ggplot(clade_data_filtered, aes(x = Mean_Geographic_Overlap, y = gamma)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_point(aes(color = Ecology, size = n.tips), alpha = 0.9) +
  geom_text_repel(aes(label = Clade_Name, color = Ecology), size = 3.5, fontface = "bold", show.legend = FALSE) +
  facet_wrap(~ Ecology, scales = "fixed") + # Fixed scales to highlight the different spatial domains
  scale_color_manual(values = c("Marine"="#1f78b4", "Terrestrial"="#33a02c", "Freshwater"="#ff7f00")) +
  labs(title = "Supp Fig 1: Overlap vs Gamma by Environment", 
       subtitle = "Highlighting the confounding effect of different dispersal capabilities (Vagility)",
       x = "Mean Geographic Overlap", y = "Gamma Statistic", size = "Number of Tips") +
  theme_bw(base_size = 14) + theme(legend.position = "none") # Legend hidden since facets are self-explanatory

plot(plot_supp_gamma)
ggsave("PLOT_3_Supplementary_Faceted_Gamma.png", plot = plot_supp_gamma, width = 10, height = 5, dpi = 300, bg = "white")

# 6B: Faceted DeltaR
plot_supp_deltar <- ggplot(clade_data_filtered, aes(x = Mean_Geographic_Overlap, y = deltaR)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_point(aes(color = Ecology, size = n.tips), alpha = 0.9) +
  geom_text_repel(aes(label = Clade_Name, color = Ecology), size = 3.5, fontface = "bold", show.legend = FALSE) +
  facet_wrap(~ Ecology, scales = "fixed") +
  scale_color_manual(values = c("Marine"="#1f78b4", "Terrestrial"="#33a02c", "Freshwater"="#ff7f00")) +
  labs(title = "Supp Fig 2: Overlap vs DeltaR by Environment", 
       subtitle = "Highlighting the confounding effect of different dispersal capabilities (Vagility)",
       x = "Mean Geographic Overlap", y = "DeltaR", size = "Number of Tips") +
  theme_bw(base_size = 14) + theme(legend.position = "none")

plot(plot_supp_deltar)
ggsave("PLOT_4_Supplementary_Faceted_DeltaR.png", plot = plot_supp_deltar, width = 10, height = 5, dpi = 300, bg = "white")

# --- 7. EXPORT FINAL DATA ----------------------------------------------------
write.csv(clade_data_filtered, "final_comparison_table.csv", row.names = FALSE)
cat("\n--- Full Analytical Pipeline Complete ---\n")
cat("Generated 5 Plots (1 Diagnostic, 2 Main, 2 Supplementary).\n")
