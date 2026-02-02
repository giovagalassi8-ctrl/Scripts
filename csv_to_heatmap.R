#!/usr/bin/env Rscript

# Description: This script is designed to visualize data from a CSV file as an annotated heatmap. It handles complex layouts by merging multi-row headers into clear column labels and applies a custom Red-White-Blue color gradient (0-50-100) to represent intensity.

# Usage: Rscript heatmap_custom.R [filename]

# IMPORTANT: PRELIMINARY CHECK OF THE .CSV FILE
# Before running this script, it is recommended to open the .csv file and carefully check the header rows (the first rows defining the columns).
# Sometimes, when converting from Excel to CSV, cells that were "merged" in Excel become one filled cell followed by empty cells in the CSV.
# If there are unwanted white/empty cells in the header rows, the script might generate incomplete labels.
# -> FILL THESE CELLS with the correct text before launching the script so that every column has a complete and correct label in the plot.


# --- 1. CONFIGURATION (CHANGE THESE VALUES AS NEEDED) ---
# DEFAULT FILE: The file to process if none is specified in terminal
DEFAULT_FILENAME <- "extracted_Unnamed_Table.csv"   

# HEADER ROWS: Number of rows at the top to use as column headers.
# Adjust the number based on the header rows present in the matrix. The specified rows will be merged to generate the labels in the chart. 
# For example, with the numbers 3 (default) the script will merge the top 3 rows into the column name.
NUM_HEADER_ROWS <- 3 

# TEXT SIZE: Size of the numbers inside the cells (0.8 is standard, lower to shrink)
TEXT_CEX <- 0.8

# COLOR SCALE LIMITS: The range of data values
MIN_VAL <- 0   # Corresponds to Red
MID_VAL <- 50  # Corresponds to White
MAX_VAL <- 100 # Corresponds to Blue

# --- 2. HELPER FUNCTIONS ---
get_user_input <- function(prompt_text) {
  if (interactive()) {
    return(readline(prompt = prompt_text))
  } else {
    cat(prompt_text)
    con <- file("stdin")
    on.exit(close(con))
    return(readLines(con, n = 1))
  }
}

prepare_data <- function(file_path, header_rows) {
  # Read entire file as text first to handle complex headers safely
  raw_lines <- readLines(file_path)
  
  # Read the data assuming the row AFTER the headers is the start of data
  # We use check.names=FALSE to prevent R from messing up duplicates temporarily
  full_df <- read.csv(text = raw_lines, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  
  # 1. Extract Header Rows
  if (header_rows > 0) {
    headers_part <- full_df[1:header_rows, , drop = FALSE]
    data_part    <- full_df[(header_rows + 1):nrow(full_df), , drop = FALSE]
    
    # Merge header rows vertically to create unique column names
    # e.g., "Model" + "MM" + "20" -> "Model.MM.20"
    # We strip empty spaces and NAs
    merged_headers <- apply(headers_part, 2, function(col) {
      clean_vals <- col[!is.na(col) & trimws(col) != ""]
      paste(clean_vals, collapse = "\n") # Use newline for cleaner plot labels
    })
    colnames(data_part) <- merged_headers
  } else {
    # No header case
    data_part <- full_df
  }
  
  # 2. Extract Row Names (Assume 1st column is the label)
  row_labels <- data_part[, 1]
  matrix_data <- data_part[, -1] # Remove label column
  
  # 3. Clean and Convert to Numeric Matrix
  # Force conversion to numeric (non-numeric becomes NA)
  mat <- suppressWarnings(sapply(matrix_data, as.numeric))
  
  # 4. Remove Empty Columns/Rows
  # Keep columns that are NOT fully NA
  valid_cols <- apply(mat, 2, function(x) !all(is.na(x)))
  mat <- mat[, valid_cols, drop = FALSE]
  
  # Keep rows that are NOT fully NA
  valid_rows <- apply(mat, 1, function(x) !all(is.na(x)))
  mat <- mat[valid_rows, , drop = FALSE]
  
  # Assign Row Names
  rownames(mat) <- row_labels[valid_rows]
  
  return(mat)
}

draw_custom_heatmap <- function(mat) {
  # Flip matrix vertically so it plots top-to-bottom (like the CSV)
  # In R image(), (0,0) is bottom-left. We want row 1 at top.
  mat_flipped <- mat[nrow(mat):1, , drop = FALSE]
  
  # Define Color Palette: Red -> White -> Blue
  # We create two gradients: Red-to-White and White-to-Blue
  n_steps <- 100
  palette_low  <- colorRampPalette(c("red", "white"))(n_steps/2)
  palette_high <- colorRampPalette(c("white", "blue"))(n_steps/2)
  full_palette <- c(palette_low, palette_high)
  
  # Calculate dimensions
  n_rows <- nrow(mat_flipped)
  n_cols <- ncol(mat_flipped)
  
  # Prepare Plot Area (Margins: Bottom, Left, Top, Right)
  # Increase left margin for row names
  par(mar = c(5, 8, 4, 2)) 
  
  # Draw Image
  # zlim ensures 0 is Red, 50 is White, 100 is Blue even if data is missing extremes
  image(1:n_cols, 1:n_rows, t(mat_flipped), 
        col = full_palette, 
        axes = FALSE, 
        xlab = "", ylab = "",
        zlim = c(MIN_VAL, MAX_VAL),
        main = "Heatmap Analysis")
  
  # Add Axes
  # X Axis (Top or Bottom? Usually Bottom for columns)
  axis(1, at = 1:n_cols, labels = colnames(mat_flipped), las = 2, cex.axis = 0.7, tick = FALSE)
  # Y Axis (Left, showing row names)
  axis(2, at = 1:n_rows, labels = rownames(mat_flipped), las = 2, cex.axis = 0.8, tick = FALSE)
  
  # Add Grid Lines (Optional, for clarity)
  grid(nx = n_cols, ny = n_rows, col = "gray90", lty = 1)
  
  # Add Values as Text
  for (x in 1:n_cols) {
    for (y in 1:n_rows) {
      val <- t(mat_flipped)[x, y]
      if (!is.na(val)) {
        text(x, y, labels = round(val, 1), col = "white", cex = TEXT_CEX, font = 2) # font=2 is bold
      }
    }
  }
}

# --- 3. MAIN EXECUTION ---

# Check arguments
args <- commandArgs(trailingOnly = TRUE)
file_path <- if (length(args) > 0) args[1] else DEFAULT_FILENAME

if (!file.exists(file_path)) stop(paste("Error: File not found:", file_path))

cat(paste("Processing file:", file_path, "with", NUM_HEADER_ROWS, "header rows...\n"))

# Load Data
final_matrix <- prepare_data(file_path, NUM_HEADER_ROWS)

# Display Stats
cat(sprintf("Matrix loaded: %d rows x %d columns.\n", nrow(final_matrix), ncol(final_matrix)))

# Draw to Screen
if (interactive() || .Platform$OS.type == "windows" || Sys.getenv("DISPLAY") != "") {
  cat("Displaying Heatmap...\n")
  # We use dev.new() to try opening a window if not in RStudio
  if (!interactive() && .Platform$OS.type == "windows") grDevices::windows()
  if (!interactive() && Sys.getenv("DISPLAY") != "") grDevices::x11()
  
  draw_custom_heatmap(final_matrix)
}

# Save Prompt
save_choice <- get_user_input("\nDo you want to save this heatmap? [y/N]: ")

if (tolower(substr(save_choice, 1, 1)) == "y") {
  clean_name <- tools::file_path_sans_ext(basename(file_path))
  out_filename <- paste0("heatmap_", clean_name, ".pdf")
  
  # Save as PDF
  pdf(out_filename, width = 12, height = 8)
  draw_custom_heatmap(final_matrix)
  dev.off()
  
  cat(paste("Saved to:", out_filename, "\n"))
} else {
  cat("Exiting.\n")
}
