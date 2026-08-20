# plot summary data from matrix stats "all_stats.csv")

library(GGally)
library(ggplot2)
library(dplyr)
library(svglite)

data <- read.csv("02_Matrix_Stats/All_matrices_stats.csv", header = TRUE, stringsAsFactors = FALSE)
data_filtered <- data %>%
  filter(grepl("_PM$", matrix))

cols_to_plot <- c("alignment.lenght", "parsimony.informative.sites", 
                  "Proportion.variable.sites", "missing.percent")

my_colors <- c("allgenes" = "skyblue4", "lb" = "grey", "rcv" = "indianred3")

lowerfun <- function(data, mapping){
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_smooth(method = "lm", se = FALSE, size = 0.8) +
    theme_bw()
}
diagfun <- function(data, mapping){
  ggplot(data = data, mapping = mapping) +
    geom_boxplot(alpha = 0.6) +
    theme_bw()
}

p <- ggpairs(data_filtered, 
             columns = cols_to_plot,
             mapping = aes(color = category, fill = category),
             lower = list(continuous = lowerfun),
             diag = list(continuous = diagfun),
             legend = 1) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(axis.text = element_text(size = 7),
        strip.text = element_text(size = 8, face = "bold"),
        legend.position = "bottom")

# save
ggsave("phylogenomic_plot.pdf", plot = p, width = 10, height = 10, device = "pdf")
ggsave("phylogenomic_plot.svg", plot = p, width = 10, height = 10, device = svglite::svglite)
