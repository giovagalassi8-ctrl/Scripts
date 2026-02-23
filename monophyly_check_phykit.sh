#!/bin/bash

# This script automates the execution of phykit monophyly_check across multiple phylogenetic tree files organized in a nested folder structure.
# It processes all .treefile files (excluding *_with_clade.treefile) found in subdirectories and consolidates the results into a single output file. 

# IMPORTANT: phykit has problems calculating parameters for groups where 2 species are present. 
#            Although it does not give an error during the process, it fails to give a result in the output.

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
# [bash] ./run_monophyly_scan.sh <group_file.txt>


# --- ARGUMENT CHECKING ---
# We ensure the user provided exactly one argument (the group file).
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <group_file.txt>"
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

# --- INITIALIZATION ---
# Create (or overwrite) the output file to start with a clean slate.
> "$OUTPUT_FILE"

echo "Starting analysis..."
echo "  > Group File: $GROUP_FILE_PATH"
echo "  > Output File: $OUTPUT_FILE"
echo "------------------------------------------------------------------"

# --- MAIN LOOP ---
# Iterate through main directories (change the folder number range if necessary - in this case the range is from 01_ to 07_).
for main_dir in 0[1-7]_*; do
    
    # Skip if it's not a directory
    [ -d "$main_dir" ] || continue

    # Iterate through sub-directories (01_ to 09_). Change folder number if necessary
    for sub_dir in "$main_dir"/0[1-9]_*; do
        
        # Skip if it's not a directory
        [ -d "$sub_dir" ] || continue

        # --- FIND THE CORRECT TREEFILE ---
        # We look for a file ending in .treefile BUT NOT _with_clade.treefile (change the -name section in the script if necessary)
        # 'find' is used with '! -name' to exclude the annotated file.
        tree_file=$(find "$sub_dir" -maxdepth 1 -type f -name "*.treefile" ! -name "*_clade.treefile" | head -n 1)      # 'head -n 1' ensures we only pick one file if multiple matches exist.

        # Check if a valid tree file was found in this folder
        if [ -n "$tree_file" ]; then
            
            # --- PREPARE OUTPUT FORMAT ---
            # Extract the filename
            base_tree_name=$(basename "$tree_file")
            
            # Remove the extension to get the clean name
            clean_name="${base_tree_name%.treefile}"

            # Write the clean name to the output file
            echo "$clean_name" >> "$OUTPUT_FILE"

            # --- DYNAMIC INTERSECTION ---
            # Define a unique temp filename for this iteration
            TEMP_LIST="temp_subset_${clean_name}.txt"
            
            # 1. sed 's/\r$//': Sanitize input file (remove Windows carriage returns)
            # 2. grep -F -o -f ... : 
            #    -F: Use fixed strings (not regex) from the group file.
            #    -f: Read patterns from the sanitized group file.
            #    -o: Output ONLY the matching part found in the tree file.
            # This effectively extracts the species names that are present in the tree.
            
            grep -F -o -f <(sed 's/\r$//' "$GROUP_FILE_PATH") "$tree_file" | sort | uniq > "$TEMP_LIST"
            
            # Count how many valid taxa were found
            MATCH_COUNT=$(wc -l < "$TEMP_LIST")
            
            # --- RUN PHYKIT ---
            # At least 2 taxa are necessary to check for monophyly
            if [ "$MATCH_COUNT" -ge 2 ]; then
                
                # Execute Phykit using the TEMPORARY list (subset), not the master list.
                # 2> /dev/null hides stderr noise.
                phykit monophyly_check "$tree_file" "$TEMP_LIST" 2> /dev/null >> "$OUTPUT_FILE"
                
                echo "[PROCESSED] $clean_name (Used $MATCH_COUNT taxa)"
                
            elif [ "$MATCH_COUNT" -eq 1 ]; then
                # Only 1 species found -> Cannot calculate monophyly, but not an error.
                echo "Single_Taxon_Present" >> "$OUTPUT_FILE"
                echo "[SKIPPED]   $clean_name (Only 1 taxon found)"
                
            else
                # 0 species found
                echo "No_Taxa_Found" >> "$OUTPUT_FILE"
                echo "[SKIPPED]   $clean_name (0 taxa found)"
            fi
            
            # CLEANUP: Remove the temporary file for this loop
            rm "$TEMP_LIST" 2>/dev/null

        else
            echo "[WARNING] No valid .treefile found in $sub_dir"
        fi
    done
done

echo "------------------------------------------------------------------"
echo "Done. Results saved in: $OUTPUT_FILE"
