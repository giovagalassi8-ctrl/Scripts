############################################################
# This script reads a generic tab-delimited text file containing summary statistics.
#   
# Object names are parsed assuming the following pattern:
#     PREFIX_A_PREFIX_B_SUFFIX.out
# Example:
#     concatenated_MS90_g80_allgenes.out
#
#   From this pattern the script extracts:
#     - Group (Y axis): MS90_g80
#     - Color group   : allgenes
#
#   Objects sharing the same group (MS90_g80) are plotted
#   on the same Y-axis line with different colors.
#
# USAGE:
#   RStudio:
#     source("generic_interactive_plot.R")
#
#   Terminal:
#     Rscript generic_interactive_plot.R
#
# REQUIRED PACKAGES:
#   install.packages(c("tidyverse"))
############################################################

library(tidyverse)

# --------------------------------------------------
# 1) INPUT FILE
# --------------------------------------------------

input_file <- "MS90_all_matrix_stats.txt"

if (!file.exists(input_file)) {
  stop(paste("ERROR: File", input_file, "not found."))
}

# --------------------------------------------------
# 2) READ DATA
# --------------------------------------------------

data_raw <- read_table(input_file)

# Rename first column as Object
colnames(data_raw)[1] <- "Object"

# --------------------------------------------------
# 3) PARSE OBJECT NAMES
# --------------------------------------------------
# Expected structure (example):
#   concatenated_MS90_g80_allgenes.out

data_parsed <- data_raw %>%
  mutate(
    Object_clean = gsub("\\.out$", "", Object),       # Removes the file extension, in this case .out; if it is different, change it accordingly.
   
    # Extract group for Y axis: after first "_" up to third "_"
    Group = sub("^[^_]+_([^_]+_[^_]+)_.*", "\\1", Object_clean),
    
    # Extract suffix for color: after third "_" to end
    ColorGroup = sub(".*_[^_]+_[^_]+_([^_]+)$", "\\1", Object_clean)
  )

# --------------------------------------------------
# 4) INTERACTIVE COLUMN SELECTION
# --------------------------------------------------

numeric_cols <- names(data_parsed)[sapply(data_parsed, is.numeric)]

cat("\nAvailable numeric columns:\n")
for (i in seq_along(numeric_cols)) {
  cat(sprintf("  [%d] %s\n", i, numeric_cols[i]))
}

cat("\nSelect the column number to plot: ")
col_choice <- as.integer(readLines(con = stdin(), n = 1))

if (is.na(col_choice) || col_choice < 1 || col_choice > length(numeric_cols)) {
  stop("Invalid selection.")
}

selected_column <- numeric_cols[col_choice]

cat("\nSelected column:", selected_column, "\n")

# --------------------------------------------------
# 5) PREPARE DATA FOR PLOTTING
# --------------------------------------------------

# Data is being processed for visualization purposes.

plot_data <- data_parsed %>%
  select(Group, ColorGroup, Value = all_of(selected_column))

# --------------------------------------------------
# 6) CREATE PLOT
# --------------------------------------------------

plot_main <- ggplot(
  plot_data,
  aes(
    x = Value,   # Numerical values will be placed on the x-axis
    y = Group,   # Groups will be placed on the y-axis
    color = ColorGroup     # Points will be colored based on the object type.
  )
) +
  geom_point(
    size = 4,    # Change this value to modify the point size in the plot.
    shape = 16   # Change this value to modify the point shape in the plot.
  ) +
  scale_x_continuous(           # Forces the x-axis to always start from zero.
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  )                                
 
 theme_minimal(
  base_size = 11,
  base_family = "",
  header_family = NULL,
  base_line_size = base_size/22,
  base_rect_size = base_size/22,
  ink = "black",
  paper = "white",
  accent = "#3366FF"
)

  labs(                      # Sets the axis labels.
    title = paste("Comparison based on", selected_column),
    x = selected_column,
    y = "Group",
    color = "Type"
  )

print(plot_main)

# --------------------------------------------------
# 7) ASK WHETHER TO SAVE
# --------------------------------------------------

cat("\nDo you want to save the plot? (y/n): ")
save_choice <- tolower(readLines(con = stdin(), n = 1))

if (save_choice == "y") {
  output_name <- paste0("plot_", selected_column, ".pdf")
  ggsave(
    filename = output_name,
    plot = plot_main,
    width = 9,
    height = 5
  )
  cat("Plot saved as:", output_name, "\n")
} else {
  cat("Plot not saved.\n")
}
