#!/bin/bash

# DESCRIPTION: This script automates the execution of phykit monophyly_check across multiple phylogenetic tree files organized in a nested folder structure. It processes all .treefile files (excluding *_with_clade.treefile) found in subdirectories and consolidates the results into a single output file. 
#
# OUTPUT: The results are aggregated into a single output file with the following format:
#         Tree_Name_Without_Extension
#         <phykit_output_result>
#
# Folder Structure Expected:
#   01_/
#     01_/*.treefile
#     02_/*.treefile
#     ...
#   02_/
#     01_/*.treefile
#     ...
#   ...
#   07_/
#
# USAGE:
#   ./run_monophyly_scan.sh <group_file.txt>


# --- 1. ARGUMENT CHECKING ---
# We ensure the user provided exactly one argument (the group file).
if [ "$#" -ne 1 ]; then
    echo "------------------------------------------------------------------"
    echo "ERROR: Missing arguments."
    echo "Usage: $0 <group_file.txt>"
    echo "Example: $0 Mollusca.txt"
    echo "------------------------------------------------------------------"
    exit 1
fi

# Store the input taxa file
GROUP_FILE_PATH="$1"

# It is extracted just the filename (e.g., "Mollusca.txt") to name the output file.
GROUP_FILENAME=$(basename "$GROUP_FILE_PATH")

# Define the output filename as requested
OUTPUT_FILE="monophyly_check_${GROUP_FILENAME}"

# Check if the provided group file actually exists
if [ ! -f "$GROUP_FILE_PATH" ]; then
    echo "ERROR: The file '$GROUP_FILE_PATH' does not exist."
    exit 1
fi

# --- 2. INITIALIZATION ---
# Create (or overwrite) the output file to start with a clean slate.
> "$OUTPUT_FILE"

echo "Starting analysis..."
echo "  > Group File: $GROUP_FILE_PATH"
echo "  > Output File: $OUTPUT_FILE"
echo "------------------------------------------------------------------"

# --- 3. MAIN LOOP ---
# Iterate through main directories (01_ to 07_). Change folder number if necessary.
for main_dir in 0[1-7]_*; do
    
    # Skip if it's not a directory
    [ -d "$main_dir" ] || continue

    # Iterate through sub-directories (01_ to 09_). Change folder number if necessary
    for sub_dir in "$main_dir"/0[1-9]_*; do
        
        # Skip if it's not a directory
        [ -d "$sub_dir" ] || continue

        # --- 4. FIND THE CORRECT TREEFILE ---
        # We look for a file ending in .treefile BUT NOT _with_clade.treefile
        # 'find' is used with '! -name' to exclude the annotated file.
        tree_file=$(find "$sub_dir" -maxdepth 1 -type f -name "*.treefile" ! -name "*_with_clade.treefile" | head -n 1)      # 'head -n 1' ensures we only pick one file if multiple matches exist.

        # Check if a valid tree file was found in this folder
        if [ -n "$tree_file" ]; then
            
            # --- 5. PREPARE OUTPUT FORMAT ---
            # Extract the filename
            base_tree_name=$(basename "$tree_file")
            
            # Remove the extension to get the clean name
            clean_name="${base_tree_name%.treefile}"

            # Write the clean name to the output file
            echo "$clean_name" >> "$OUTPUT_FILE"

            # --- 6. RUN PHYKIT ---
            # Execute the analysis and append the result to the next line of the output file
            phykit monophyly_check "$tree_file" "$GROUP_FILE_PATH" >> "$OUTPUT_FILE"
            
            # Print a small progress indicator to the screen
            echo "[PROCESSED] $clean_name"

        else
            # Warning if no suitable tree file is found in a subdirectory
            echo "[WARNING] No valid .treefile found in $sub_dir"
        fi

    done
done

echo "------------------------------------------------------------------"
echo "Analysis complete. Results saved in: $OUTPUT_FILE"
