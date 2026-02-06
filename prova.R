#!/usr/bin/env Rscript

# --- CONFIGURATION ---
# Input file
INPUT_FILE  <- "extracted_Unnamed_Table.csv"

# Number of rows to merge for the header (e.g., 3 rows become one label)
HEADER_ROWS <- 3 

# VISUAL SETTINGS
# Increase GRID_WIDTH to make squares smaller and margins larger (Try 3, 5, or 10)
GRID_WIDTH  <- 5  
TEXT_SIZE   <- 0.7

# Color scale (Red=0, White=50, Blue=100)
MIN_VAL <- 0
MID_VAL <- 50
MAX_VAL <- 100

# Output filename (PDF)
OUTPUT_FILE <- "heatmap_output.pdf"

# --- MAIN SCRIPT ---

# 1. Read Data
if (!file.exists(INPUT_FILE)) {
  stop(paste("Error: File not found:", INPUT_FILE))
}

# Read raw lines to handle irregular headers
raw_lines <- readLines(INPUT_FILE)
full_df   <- read.csv(text = raw_lines, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)

# 2. Process Headers
if (HEADER_ROWS > 0) {
  # Extract top rows for headers
  headers_part <- full_df[1:HEADER_ROWS, , drop = FALSE]
  data_part    <- full_df[(HEADER_ROWS + 1):nrow(full_df), , drop = FALSE]
  
  # Merge header columns vertically
  merged_headers <- apply(headers_part, 2, function(col) {
    # Remove NAs and empty strings, then join with a newline
    clean_vals <- col[!is.na(col) & trimws(col) != ""]
    paste(clean_vals, collapse = "\n") 
  })
  colnames(data_part) <- merged_headers
} else {
  data_part <- full_df
}

# 3. Process Data Matrix
# Assume first column is Row Names
row_labels <- data_part[, 1]
matrix_data <- data_part[, -1]

# Convert to numeric matrix (suppress warnings for NAs)
mat <- suppressWarnings(sapply(matrix_data, as.numeric))
rownames(mat) <- row_labels

# Remove rows/cols that are completely empty/NA
mat <- mat[, colSums(!is.na(mat)) > 0, drop = FALSE]
mat <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]

# Flip matrix (so row 1 is at the top in the plot)
mat_flipped <- mat[nrow(mat):1, , drop = FALSE]

# 4. Define Colors (Red -> White -> Blue)
steps <- 100
palette <- c(
  colorRampPalette(c("red", "white"))(steps/2),
  colorRampPalette(c("white", "blue"))(steps/2)
)

# 5. Plotting Function
draw_plot <- function() {
  n_rows <- nrow(mat_flipped)
  n_cols <- ncol(mat_flipped)
  
  # Margins: Bottom, Left, Top, Right
  par(mar = c(5, 8, 4, 2))
  
  # Draw base heatmap
  image(1:n_cols, 1:n_rows, t(mat_flipped), 
        col = palette, 
        axes = FALSE, 
        xlab = "", ylab = "",
        zlim = c(MIN_VAL, MAX_VAL),
        main = "Heatmap Analysis")
  
  # --- THE TRICK: Draw thick white grid lines to create separation ---
  # Vertical white lines
  abline(v = seq(0.5, n_cols + 0.5, 1), col = "white", lwd = GRID_WIDTH)
  # Horizontal white lines
  abline(h = seq(0.5, n_rows + 0.5, 1), col = "white", lwd = GRID_WIDTH)
  # -----------------------------------------------------------------
  
  # Add Labels
  axis(1, at = 1:n_cols, labels = colnames(mat_flipped), las = 2, tick = FALSE, cex.axis = 0.7)
  axis(2, at = 1:n_rows, labels = rownames(mat_flipped), las = 2, tick = FALSE, cex.axis = 0.8)
  
  # Add Numbers (plotted ON TOP of the white grid)
  for (x in 1:n_cols) {
    for (y in 1:n_rows) {
      val <- t(mat_flipped)[x, y]
      if (!is.na(val)) {
        # Check if value is dark (blue/red) to use white text, else black
        text_col <- ifelse(abs(val - MID_VAL) > 25, "white", "black")
        text(x, y, labels = round(val, 1), col = text_col, cex = TEXT_SIZE, font = 2)
      }
    }
  }
}

# 6. Execute and Save
cat("Generating PDF...\n")
pdf(OUTPUT_FILE, width = 12, height = 8)
draw_plot()
dev.off()

cat(paste("Done! Saved to", OUTPUT_FILE, "\n"))
