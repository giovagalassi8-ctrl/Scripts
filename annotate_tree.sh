#!/bin/bash

# This script annotates a phylogenetic tree file by appending a clade suffix to the species names, based on a provided mapping table.

# USAGE: ./annotate_tree.sh <mapping_table.tsv> <tree_file>

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Error: Invalid number of arguments."
    echo "Usage: $0 <mapping_table.tsv> <tree_file>"
    exit 1
fi

# Assign arguments to readable variables
MAPPING_FILE="$1"      # Mapping file (tsv format) contains species names in column 1 and their corresponding clades in column 2.
TREE_FILE="$2"         # The target tree file where labels will be added.

# Generate output filename
# We take the full input filename ($TREE_FILE) and append the suffix directly.
BASE_NAME="${TREE_FILE%.*}"
OUTPUT_FILE="${BASE_NAME}_annotated.treefile"

echo "Starting annotation process..."
echo "  > Reading mapping from: $MAPPING_FILE"
echo "  > Reading tree from:    $TREE_FILE"

# Start awk processing
awk -v OFS="" '
    # --- PHASE 1: Load Mapping Data ---
    NR==FNR {
        # Extracts the first 3 letters of the Clade name to create a suffix.
        clade3[$1] = substr($2, 1, 3)
        next 
    }

    # --- PHASE 2: Process Tree File ---
    # Scans the provided tree file and replaces every known species name with "SpeciesName_Suffix".
    {
        for (species in clade3) {
            gsub(species, species"_"clade3[species])
        }
        print
    }
' "$MAPPING_FILE" "$TREE_FILE" > "$OUTPUT_FILE"

echo "Done! The annotated tree has been saved as: $OUTPUT_FILE"
