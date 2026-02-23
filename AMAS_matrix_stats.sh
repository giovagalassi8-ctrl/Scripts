#!/bin/bash

#This script recursively finds all concatenated_*.out files, copies them into the 11_matrix_stats directory, and runs AMAS.py to generate summary statistics for each file, producing one summary output per concatenated matrix.

# USAGE:
# [bash] ./AMAS_matrix_stats.sh

# Output directory where copied files and AMAS results will be stored
OUTDIR="11_matrix_stats"

# Create output directory if it does not exist
mkdir -p "$OUTDIR"

# Recursively find all concatenated_*.out files
find . -type f -name "concatenated_*.out" | while read -r file; do
    # Extract filename without path
    basefile=$(basename "$file")

    # Extract the part between 'concatenated_' and '.out'
    core_name="${basefile#concatenated_}"
    core_name="${core_name%.out}"

    # Destination path for the copied concatenated file
    copied_file="$OUTDIR/$basefile"

    # Copy the file (overwrite if it already exists)
    cp "$file" "$copied_file"

    # Name of the AMAS summary output file
    summary_file="$OUTDIR/summary_${core_name}.txt"

    # Run AMAS summary on the concatenated file (add the correct AMAS.py path if necessary)
    python3 AMAS.py summary \
        -i "$copied_file" \
        -f fasta \
        -d aa \
        -o "$summary_file"

    echo "Processed: $basefile → summary_${core_name}.txt"
done
