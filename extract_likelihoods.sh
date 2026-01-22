#!/bin/bash

# Define the output filename
output_file="Likelyhood.txt"

# Initialize/Clear the output file
> "$output_file"

echo "Starting extraction. Looking for Base_results.txt or Gamma_results.txt..."

# Find files matching either name.
# We use \( ... -o ... \) to group the conditions (OR logic).
find . -type f \( -name "Base_results.txt" -o -name "Gamma_results.txt" \) | sort | while read -r filepath; do

    # Extract the directory path (e.g., ./00_1L/1K/1N)
    dir_path=$(dirname "$filepath")

    # Clean the path by removing the leading "./"
    clean_path=${dir_path#./}

    # Extract the first line of the found file
    first_line=$(head -n 1 "$filepath")

    # Write the directory path and the extracted line to the output file
    # We use a tab (\t) as separator
    echo -e "${clean_path}\t${first_line}" >> "$output_file"

done

echo "Process completed. Data saved to $output_file"
