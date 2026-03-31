# !/usr/env/Rscript

# This script allows you to create a scatterplot to graphically represent the results of 'phykit toverr'.
# In particular, the treeness (measure of the proportion of the total tree length (sum of all branch lengths) that is found on internal branches) is represented on the x-axis, 
# while the Relative Composition Variability (RCV, measure of the average variability in sequence composition among taxa in a sequence alignment) is represented on the y-axis.
# Higher values of treeness/RCV (high treeness, low RCV) are desirable, as they indicate that a gene is likely to be less susceptible to systematic biases.




library(ggplot2)

data <- read.table('05_Treeness/allgenes_treeness.tsv', header = FALSE)
data <- data[,-2]

colnames(data)[2] <- "treeness"
colnames(data)[3] <- "RCV"

kpic <- data[c(16:18),]

col <- ifelse(data$V1 %in% gappyout$V1, "red", "lightgrey")

p <- ggplot(data, aes(
    x=treeness, 
    y=RCV)
    )+
  geom_point(
    size = 2.5,
    alpha = 0.85,
    col = col
    )+
  geom_abline(
    intercept = 0,
    slope = 1,
    )+
  coord_fixed(ratio = 1)+
  theme_bw()+
  scale_x_continuous(limits = c(0.215, 0.255))+
  scale_y_continuous(limits = c(0.215, 0.255))+
  labs(
    x = "Treeness",
    y = "Relative Composition Variability (RCV)",
    col = NULL
  )+
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    axis.text = element_text(color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 12),
  )
