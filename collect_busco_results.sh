#!/bin/bash

# This script automates the aggregation of BUSCO results from multiple directories. 
# It scans the current folder for BUSCO output directories, identifies the specific lineage run (e.g., metazoa, arthropoda), and copies the result files into a single consolidated directory.

# USAGE: 
# [bash] ./collect_busco_results.sh
# (Run it in the directory containing your BUSCO output folders)

# --- CONFIGURATION VARIABLES ---
# To change the input suffix or output folder, edit the "CONFIGURATION" section below before running.

# Looks for directories matching a specific suffix (CHANGE THE SUFFIX IF NECESSARY; default: "_busco")
INPUT_SUFFIX="_busco"

# The name of the new directory where files will be collected (CHANGE THE NAME OF THE RESULTING DIRECTORY IF NECESSARY)
OUTPUT_DIR="all_busco_sequences"

# --- MAIN SCRIPT ---

# Create the main output directory if it doesn't exist
# -p ensures no error if it already exists
mkdir -p "$OUTPUT_DIR"

echo "Starting BUSCO collection..."
echo "----------------------------------------"

# Iterate through all directories ending with the specified suffix
for dir in *"$INPUT_SUFFIX"; do

    # Check if the match is actually a directory (handles cases where no files match)
    if [ -d "$dir" ]; then
        
        # 1. Extract the Species/Sample Name
        # We use bash parameter expansion (%) to remove the suffix (faster and cleaner than piping to sed).
        species_name="${dir%$INPUT_SUFFIX}"
        
        # 2. Dynamically find the run directory (Lineage agnostic)
        # BUSCO creates a folder named 'run_[lineage_name]'.
        # We search for any item starting with 'run_' inside the current busco dir.
        # We assume there is only one run folder per busco directory.
        run_folder=$(find "$dir" -maxdepth 1 -type d -name "run_*" | head -n 1)

        # Check if a run folder was actually found and contains the sequence subdir
        if [ -n "$run_folder" ] && [ -d "$run_folder/busco_sequences" ]; then
            
            echo "Processing: $species_name"
            
            # Create the species-specific destination folder
            target_path="$OUTPUT_DIR/$species_name"
            mkdir -p "$target_path"

            # Copy the run folder (recursively) to the destination
            # We use cp -r to copy the directory structure
            cp -r "$run_folder" "$target_path/"
            
        else
            # Print a warning if the expected structure is missing
            # >&2 redirects echoes to standard error (good practice for logs)
            echo "WARNING: No valid 'busco_sequences' found in $dir" >&2
        fi
        
    else
        # If the loop doesn't find any directory matching the pattern
        echo "No directories ending in '$INPUT_SUFFIX' found in $(pwd)."
        break
    fi

done

echo "----------------------------------------"
echo "Job finished. All files are in: $OUTPUT_DIR"
