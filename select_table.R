#!/usr/bin/env Rscript

# Description: Parses a CSV file containing multiple tables. It detects tables stacked vertically (separated by empty lines) AND tables placedside-by-side (separated by empty columns). Allows interactive selection and saving.

# Usage:       
# [bash] Rscript select_table.R [filename]
# [RStudio] source("select_table.R")


# --- Configuration ---
# Default filename to look for if no argument is provided
DEFAULT_FILENAME <- "02_SPIRALIA_DATASET.xlsx - matrices.csv"     # Change file name here

# --- Functions ---
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

# Checks if a vector is "empty" (all NAs or empty strings)
is_column_empty <- function(x) {
  # Treat NA and empty strings as empty
  x <- as.character(x)
  all(is.na(x) | trimws(x) == "")
}

extract_subtables <- function(lines, block_start, block_end) {
  # 1. Create a raw dataframe from the lines to analyze columns
  # We use header=FALSE to treat the first row as data for detecting empty cols
  txt_conn <- textConnection(lines[block_start:block_end])
  raw_df <- tryCatch({
    read.csv(txt_conn, header = FALSE, stringsAsFactors = FALSE, colClasses = "character")
  }, error = function(e) return(NULL))
  close(txt_conn)
  
  if (is.null(raw_df) || ncol(raw_df) == 0) return(list())
  
  # 2. Identify empty columns
  # Returns a boolean vector: TRUE if column is empty, FALSE contains data
  empty_cols_mask <- sapply(raw_df, is_column_empty)
  
  # 3. Identify ranges of consecutive NON-empty columns
  # We use Run Length Encoding (rle) on the inverse mask
  has_data <- !empty_cols_mask
  rle_res <- rle(has_data)
  
  subtables <- list()
  current_col_idx <- 1
  
  for (i in seq_along(rle_res$lengths)) {
    is_data_chunk <- rle_res$values[i]
    chunk_len <- rle_res$lengths[i]
    
    if (is_data_chunk) {
      # This is a valid table block
      col_start <- current_col_idx
      col_end <- current_col_idx + chunk_len - 1
      
      # Extract the subset
      sub_df <- raw_df[, col_start:col_end, drop = FALSE]
      
      # Try to extract a name from the top-left cell
      name_candidate <- sub_df[1, 1]
      if (is.na(name_candidate) || name_candidate == "") name_candidate <- "Unnamed Table"
      
      subtables[[length(subtables) + 1]] <- list(
        name = name_candidate,
        data = sub_df,
        orig_rows = paste0(block_start, "-", block_end),
        orig_cols = paste0(col_start, "-", col_end)
      )
    }
    
    current_col_idx <- current_col_idx + chunk_len
  }
  
  return(subtables)
}

parse_all_tables <- function(lines) {
  # 1. Identify Vertical Blocks (empty lines)
  is_separator <- grepl("^,*\\s*$", lines)
  
  vertical_blocks <- list()
  current_start <- NA
  
  for (i in seq_along(lines)) {
    if (!is_separator[i]) {
      if (is.na(current_start)) current_start <- i
    } else {
      if (!is.na(current_start)) {
        vertical_blocks[[length(vertical_blocks) + 1]] <- c(start = current_start, end = i - 1)
        current_start <- NA
      }
    }
  }
  if (!is.na(current_start)) {
    vertical_blocks[[length(vertical_blocks) + 1]] <- c(start = current_start, end = length(lines))
  }
  
  # 2. Process each vertical block to find Horizontal splits
  all_tables <- list()
  
  for (v_block in vertical_blocks) {
    found_tables <- extract_subtables(lines, v_block['start'], v_block['end'])
    all_tables <- c(all_tables, found_tables)
  }
  
  return(all_tables)
}

# --- Main Execution ---
args <- commandArgs(trailingOnly = TRUE)
file_path <- if (length(args) > 0) args[1] else DEFAULT_FILENAME

if (!file.exists(file_path)) stop(paste("Error: File not found:", file_path))

cat(paste("Reading file:", file_path, "...\n"))
all_lines <- readLines(file_path, warn = FALSE)

# Detect all tables (vertical + horizontal)
tables <- parse_all_tables(all_lines)

if (length(tables) == 0) {
  stop("No tables found in the file.")
}

# Menu Selection
selected_index <- 1

if (length(tables) > 1) {
  cat("\n--- Tables Detected ---\n")
  for (i in seq_along(tables)) {
    tbl <- tables[[i]]
    cat(sprintf("[%d] %s (Rows: %s, Cols: %s)\n", 
                i, tbl$name, tbl$orig_rows, tbl$orig_cols))
  }
  
  repeat {
    input <- get_user_input("\nEnter table number: ")
    choice <- as.integer(input)
    
    if (!is.na(choice) && choice >= 1 && choice <= length(tables)) {
      selected_index <- choice
      break
    } else {
      cat("Invalid selection.\n")
    }
  }
} else {
  cat("\nOnly one table found. Auto-selecting.\n")
}

# Get selected table
sel_table_obj <- tables[[selected_index]]
final_df <- sel_table_obj$data

cat(paste("\nSelected:", sel_table_obj$name, "\n"))

# --- Header Adjustment ---
# Since we read with header=FALSE to detect blocks, the first row is now row 1 of data.
# Often in these stacked files, the first row is the 'Title' (e.g. "MS90") and 
# the SECOND row is the actual header (e.g. "model, MM, LG...").
# Let's ask the user or use a heuristic. 
# Heuristic: If row 1 has only 1 non-empty cell and row 2 has many, row 2 is likely header.

cat("--- Processing Headers ---\n")
# Simple default: assume row 2 is header if row 1 looks like a title (mostly empty)
# Otherwise assume row 1 is header.
non_empty_r1 <- sum(final_df[1,] != "")
non_empty_r2 <- sum(final_df[2,] != "")

use_row2_as_header <- FALSE
if (nrow(final_df) > 1) {
  # If row 1 has significantly fewer filled cells than row 2, it's likely a title
  if (non_empty_r1 < (non_empty_r2 / 2) || non_empty_r1 == 1) {
    use_row2_as_header <- TRUE
  }
}

if (use_row2_as_header) {
  cat("Detected Title in Row 1. Using Row 2 as Header.\n")
  colnames(final_df) <- final_df[2, ]
  final_df <- final_df[-c(1, 2), ] # Remove title and header rows
} else {
  cat("Using Row 1 as Header.\n")
  colnames(final_df) <- final_df[1, ]
  final_df <- final_df[-1, ]
}

# Reset row names
rownames(final_df) <- NULL

# Preview
cat("\n--- Data Preview ---\n")
print(head(final_df, 5))

# Export Option
save_prompt <- get_user_input("\nSave to CSV? [y/N]: ")
if (tolower(substr(save_prompt, 1, 1)) == "y") {
  clean_name <- gsub("[^a-zA-Z0-9]", "_", sel_table_obj$name)
  out_name <- paste0("extracted_", clean_name, ".csv")
  write.csv(final_df, out_name, row.names = FALSE)
  cat(paste("Saved to:", out_name, "\n"))
}

invisible(final_df)
