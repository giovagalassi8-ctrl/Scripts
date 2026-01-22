#!/bin/bash

# Define the output filename as requested
output_file="Likelyhood.txt"

# Initialize the output file (clears it if it already exists)
> "$output_file"

echo "Starting extraction of first lines from Base_results.txt..."

# Find all instances of Base_results.txt in subdirectories
# We use 'sort' to ensure the output order is consistent (e.g., 00_1L before 00_2L)
find . -name "Base_results.txt" | sort | while read -r filepath; do

    # Extract the directory path relative to the current folder
    # Example: ./00_1L/1K/1N/Base_results.txt -> ./00_1L/1K/1N
    dir_path=$(dirname "$filepath")

    # Remove the leading "./" for a cleaner output format
    clean_path=${dir_path#./}

    # Extract the very first line of the file
    first_line=$(head -n 1 "$filepath")

    # Write the path and the extracted line to the output file
    # Using a tab (\t) to separate the path from the value for readability
    echo -e "${clean_path}\t${first_line}" >> "$output_file"

done

echo "Process completed. Data saved to $output_file"
