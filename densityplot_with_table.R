# This script visualizes phylogenetic tree metrics by generating combined density and dot plots for Mean Branch Length, Sum Branch Length, and Mean Support. It reads data from text files, 
# builds individual plots with custom themes, and finally combines them into a single layout. Optionally, it computes and appends summary statistics tables next to each corresponding plot.

library(ggplot2)
library(patchwork)
library(gridExtra)

# Create all the objects for the density plot.
meanbrlength <- read.table("gotree_meanbrlenght.txt", header = FALSE, col.names = c("Model", "MeanBranchLength"))
sumbrlength <- read.table("gotree_sumbrlenght.txt", header = FALSE, col.names = c("Model", "SumBranchLength"))
meansupport <- read.table("gotree_meansupport.txt", header = FALSE, col.names = c("Model", "MeanSupport"))

# Create the density plots.
density_plot <- ggplot(meanbrlength, aes(x = MeanBranchLength)) +
  geom_density( aes(y = after_stat(scaled)*0.4), fill = "#008080", alpha = 0.4, color = "darkslategray", bounds = c(0, 1)) +
  geom_dotplot(fill = "black", alpha = 0.8, binwidth = 0.005, dotsize = 1, method = "histodot") +
  geom_vline(aes(xintercept = mean(MeanBranchLength, na.rm = TRUE)), color = "white", linetype = "solid", linewidth = 0.8) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(
    title = "Mean Branch Length Distribution",
    x = "Mean Branch Length",
    y = "Density"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12)
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1))

density_plot2 <- ggplot(sumbrlength, aes(x = SumBranchLength)) +
  geom_density( aes(y = after_stat(scaled)*0.4), fill = "#008080", alpha = 0.4, color = "darkslategray", bounds = c(0, 100)) +
  geom_dotplot(fill = "black", alpha = 0.8, binwidth = 0.5, dotsize = 1, method = "histodot") +
  geom_vline(aes(xintercept = mean(SumBranchLength, na.rm = TRUE)), color = "White", linetype = "solid", linewidth = 0.8) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(
    title = "Sum Branch Length Distribution",
    x = "Sum Branch Length",
    y = "Density"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12)
  ) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10))

density_plot3 <- ggplot(meansupport, aes(x = MeanSupport)) +
  geom_density( aes(y = after_stat(scaled)*0.6), fill = "#008080", alpha = 0.4, color = "darkslategray", bounds = c(0, 100)) +
  geom_dotplot(fill = "black", alpha = 0.8, binwidth = 0.5, dotsize = 1, method = "histodot") +
  geom_vline(aes(xintercept = mean(MeanSupport, na.rm = TRUE)), color = "White", linetype = "solid", linewidth = 0.8) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(
    title = "Mean Support Distribution",
    x = "Mean Support",
    y = "Density"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12)
  ) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10))

# Merge all the plots in one.
final_plot <- density_plot / density_plot2 / density_plot3

# See the result.
final_plot



# OPTIONALLY: The next step is to add tables containing summary statistics to the side of each chart.
# Create a data frame containing summary statistics.
create_summary_table <- function(numeric_vector) {
  stats <- summary(numeric_vector)
  df <- data.frame(
    Statistic = names(stats),
    Value = round(as.numeric(stats), 3) # Round to 3 decimals for readability
  )
  return(df)
}

# Create a theme to shrink the summary table.
small_theme <- ttheme_default(base_size = 8, padding = unit(c(2,2), "mm"))

# Convert the summary data frames into graphical objects.
# Setting rows = NULL removes the default row numbers from the left side.
table1_grob <- tableGrob(create_summary_table(meanbrlength$MeanBranchLength), rows = NULL, theme = small_theme)
table2_grob <- tableGrob(create_summary_table(sumbrlength$SumBranchLength), rows = NULL, theme = small_theme)
table3_grob <- tableGrob(create_summary_table(meansupport$MeanSupport), rows = NULL, theme = small_theme)

# Combine each plot with its corresponding table horizontally.
row1 <- density_plot + table1_grob + plot_layout(widths = c(4, 1))
row2 <- density_plot2 + table2_grob + plot_layout(widths = c(4, 1))
row3 <- density_plot3 + table3_grob + plot_layout(widths = c(4, 1))

# Merge all the rows vertically.
final_layout <- row1 / row2 / row3

# See the result.
final_layout
