library(ggplot2)
library(tidyverse)
library(viridis)

data <- read.csv("00_Data/Literature - Copy.csv", header=FALSE, stringsAsFactors = FALSE)

transposed_data <- as.data.frame(t(data), stringsAsFactors = FALSE)

objects <- as.character(unlist(data[1,-1]))
time <- as.numeric(unlist(data[2,-1]))
phyla <- as.numeric(unlist(data[4,-1]))
number_of_species <- as.numeric(unlist(data[3,-1]))

data_plot <- data.frame(
  Objects = objects,
  Time = time,
  Phyla = phyla
)


lollipop_plot <- ggplot(data_plot, aes(x = Time, y = Phyla)) +
  geom_segment(aes(x = Time, xend = Time, y = 0, yend = Phyla), color = "gray50", linewidth = 1) +
  geom_point(size = 4, color = "steelblue") +
  geom_text(aes(label = Objects), vjust = -1.5, size = 3.5, color = "black") +
  labs(
    title = "Lollipop Plot",
    x = "Timeline",
    y = "Number of Phyla"
  ) +
  theme_bw() +
  scale_y_continuous(limits = c(0,30), breaks = seq(0,30, by = 5))
