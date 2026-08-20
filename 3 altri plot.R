library(ggplot2)
library(patchwork)

# pattern: ML_MS{80|90}_{trimming}_{category}_{6aa|MM|PM}

data <- read.csv("02_Matrix_Stats/All_matrices_stats.csv", stringsAsFactors = FALSE)
data$supermatrix <- ifelse(grepl("MS80", data$matrix), "MS80", "MS90")
data$approach <- NA
data$approach[grepl("_6aa$", data$matrix)] <- "Dayhoff-6"
data$approach[grepl("_MM$",  data$matrix)] <- "Mixture"
data$approach[grepl("_PM$",  data$matrix)] <- "Partition"
# as facter
data$category    <- factor(data$category,    levels = c("allgenes", "rcv", "lb"))
data$supermatrix <- factor(data$supermatrix, levels = c("MS80", "MS90"))
data$approach    <- factor(data$approach,    levels = c("Partition", "Mixture", "Dayhoff-6"))

cat_colors <- c("allgenes" = "blue",
                "rcv"      = "red",
                "lb"       = "darkgrey")
shape_map <- c("MS80" = 16, "MS90" = 17)
alpha_map <- c("Partition" = 1, "Mixture" = 1, "Dayhoff-6" = 0.5)
theme_clean <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey95", color = NA))
params <- c("alignment.lenght",
            "parsimony.informative.sites",
            "Proportion.variable.sites",
            "missing.percent",
            "Treeness",
            "RCV",
            "meansupport",
            "sumbrlen")
param_labels <- c("alignment.lenght" = "Alignment length (bp)",
                  "parsimony.informative.sites" = "Parsimony-informative sites",
                  "Proportion.variable.sites" = "Proportion variable sites",
                  "missing.percent" = "Missing data (%)",
                  "Treeness" = "Treeness",
                  "RCV" = "RCV",
                  "meansupport" = "Mean bootstrap support",
                  "sumbrlen" = "Sum of branch lengths")

plot_list <- lapply(params, function(p) {
  ggplot(data, aes(x = category, y = .data[[p]])) +
    geom_boxplot(aes(fill = category, color = category),
                 alpha = 0.15, outlier.shape = NA, width = 0.6) +
    geom_jitter(aes(color = category, shape = supermatrix, alpha = approach),
                position = position_jitter(width = 0.18, seed = 42),
                size = 2) +
    scale_fill_manual(values = cat_colors,  guide = "none") +
    scale_color_manual(values = cat_colors, guide = "none") +
    scale_shape_manual(values = shape_map) +
    scale_alpha_manual(values = alpha_map) +
    labs(x = NULL, y = param_labels[[p]]) +
    theme_clean +
    theme(legend.position = "none")
})
# same legend
legend_plot <- ggplot(data, aes(x = category, y = Treeness)) +
  geom_point(aes(color = category, shape = supermatrix, alpha = approach), size = 3) +
  scale_color_manual(values = cat_colors, name = "Filtering scheme") +
  scale_shape_manual(values = shape_map,  name = "Supermatrix") +
  scale_alpha_manual(values = alpha_map,  name = "Approach") +
  theme_clean +
  theme(legend.position = "bottom")

shared_legend <- cowplot::get_legend(legend_plot)
boxplot_grid <- wrap_plots(plot_list, ncol = 4) /
  shared_legend +
  plot_layout(heights = c(20, 1))
scatter_plot <- ggplot(data,
                       aes(x = RCV, y = Treeness,
                           color = category, shape = supermatrix, alpha = approach)) +
  geom_point(size = 3) +
  scale_color_manual(values = cat_colors, name = "Filtering scheme") +
  scale_shape_manual(values = shape_map,  name = "Supermatrix") +
  scale_alpha_manual(values = alpha_map,  name = "Approach") +
  labs(x = "RCV (relative compositional variability)",
       y = "Treeness") +
  theme_clean +
  theme(legend.position = "bottom")

outdir <- "./plots"
dir.create(outdir, showWarnings = FALSE)

for (ext in c("pdf", "svg", "png")) {
  ggsave(file.path(outdir, paste0("matrix_boxplots.", ext)),
         boxplot_grid, width = 13, height = 7, dpi = 300)
  ggsave(file.path(outdir, paste0("treeness_vs_rcv.", ext)),
         scatter_plot, width = 13, height = 7, dpi = 300)
}
