#!/usr/bin/env Rscript

# Description: Parses a CSV file containing multiple stacked tables (separated by empty lines), allows the user to interactively select one, and loads/saves the selected table.

# Usage:       
# Rscript select_table.R [filename]
# In RStudio: source("select_table.R")


# --- Configuration ---
# Default filename to look for if no argument is provided
DEFAULT_FILENAME <- "02_SPIRALIA_DATASET.xlsx - matrices.csv"     # Change file name here

# --- Functions ---

get_user_input <- function(prompt_text) {
  # Handles input reading for both interactive session (RStudio) and non-interactive (Terminal)
  if (interactive()) {
    return(readline(prompt = prompt_text))
  } else {
    cat(prompt_text)
    # Read from stdin
    con <- file("stdin")
    on.exit(close(con))
    return(readLines(con, n = 1))
  }
}

parse_blocks <- function(lines) {
  # Identifies start and end indices of data blocks separated by empty lines
  # Returns a list of lists with 'start', 'end', and 'name'
  
  # Regex for empty or comma-only lines
  is_separator <- grepl("^,*\\s*$", lines)
  
  blocks <- list()
  current_start <- NA
  
  for (i in seq_along(lines)) {
    if (!is_separator[i]) {
      if (is.na(current_start)) {
        current_start <- i
      }
    } else {
      if (!is.na(current_start)) {
        blocks[[length(blocks) + 1]] <- c(start = current_start, end = i - 1)
        current_start <- NA
      }
    }
  }
  # Capture the last block if file doesn't end with a newline
  if (!is.na(current_start)) {
    blocks[[length(blocks) + 1]] <- c(start = current_start, end = length(lines))
  }
  
  return(blocks)
}

extract_block_name <- function(first_line) {
  # Extracts a potential name from the first cell of the block
  parts <- strsplit(first_line, ",")[[1]]
  # Clean up quotes
  parts <- gsub('"', '', parts)
  # Find first non-empty part
  name <- parts[parts != ""][1]
  if (is.na(name) || name == "") return("Unnamed Table")
  return(name)
}

# --- Main Execution ---

# Determine file path
args <- commandArgs(trailingOnly = TRUE)
file_path <- if (length(args) > 0) args[1] else DEFAULT_FILENAME

# Check if file exists
if (!file.exists(file_path)) {
  stop(paste("Error: File not found:", file_path))
}

cat(paste("Reading file:", file_path, "...\n"))
all_lines <- readLines(file_path, warn = FALSE)

# Parse blocks
blocks <- parse_blocks(all_lines)

if (length(blocks) == 0) {
  stop("No tables found in the file.")
}

# Extract names for menu
block_names <- sapply(blocks, function(b) extract_block_name(all_lines[b['start']]))

# Select table
selected_index <- 1

if (length(blocks) > 1) {
  cat("\n--- Multiple Tables Found ---\n")
  for (i in seq_along(block_names)) {
    cat(sprintf("[%d] %s (Rows %d-%d)\n", i, block_names[i], blocks[[i]]['start'], blocks[[i]]['end']))
  }
  
  repeat {
    input <- get_user_input("\nEnter the number of the table to process: ")
    choice <- as.integer(input)
    
    if (!is.na(choice) && choice >= 1 && choice <= length(blocks)) {
      selected_index <- choice
      break
    } else {
      cat("Invalid selection. Please try again.\n")
    }
  }
} else {
  cat("\nOnly one table found. Selecting automatically.\n")
}

# Process selected block
sel_block <- blocks[[selected_index]]
cat(paste("\nSelected:", block_names[selected_index], "\n"))

# Extract lines for the selected table
table_lines <- all_lines[sel_block['start']:sel_block['end']]

# Attempt to load into DataFrame
# We assume the table format allows for CSV parsing. 
# Note: Complex headers (multi-row) might require manual adjustment later.
df <- read.csv(text = table_lines, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# Show preview
cat("\n--- Data Preview (Top 5 rows) ---\n")
print(head(df, 5))
cat(sprintf("\nDimensions: %d rows x %d columns\n", nrow(df), ncol(df)))

# Option to save to file (useful for Bash workflow)
save_prompt <- get_user_input("\nDo you want to save this table to a new CSV file? [y/N]: ")
if (tolower(substr(save_prompt, 1, 1)) == "y") {
  out_name <- paste0("extracted_", gsub("[^a-zA-Z0-9]", "_", block_names[selected_index]), ".csv")
  write.csv(df, out_name, row.names = FALSE)
  cat(paste("Saved to:", out_name, "\n"))
}

# If running in RStudio/interactive mode, return the dataframe invisibly so it can be assigned
# e.g., my_data <- source("select_table.R")$value
invisible(df)
