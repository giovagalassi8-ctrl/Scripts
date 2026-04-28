# !/bin/bash

# This script performs a monophyly check using the Monophylo.py script. 
# It must be launched inside a folder containing all the treefiles you want to analyze, using a TSV file that contains the various taxonomic groupings for each species considered. 
# In this last file can be added columns also containing taxonomic groups such as superphyla, superfamilies, or any type of grouping desired

# USAGE:
# [bash] ./run_MonoPhylo.sh

# Define the output directory.
OUTDIR="OUTPUT_DIRECTORY_PATH"
# Define the absolute path for the directory containing the trees and mapping file.
CHECKDIR="DIRECTORY_PATH"
# Define the path for the final aggregated results file.
FINAL="$OUTDIR/FILE_NAME.TXT"
# Define a temporary directory for each run's intermediate outputs.
TMPDIR="$OUTDIR/tmp_run"

# Write the TSV header to the final output file (overwrites if it already exists).
echo -e "source_tree\tGrouping\tNumber_Contained_Taxa\tMonophyletic\tCategory\tSupport\tNumber_Interfering_Species\tInterfering_Species" > "$FINAL"

# Iterate over all .treefile files in the current working directory.
for tree in *.treefile; do
  # Create the temporary directory if it does not exist.
  mkdir -p "$TMPDIR"
  # Execute MonoPhylo for monophyly analysis
  python3 MonoPhylo.py \
    --tree "$CHECKDIR/$tree" \
    --out_dir "$TMPDIR" \
    --map "$CHECKDIR/TAXONOMY_FILE.CSV" \
    --support

  # Iterate through all generated text files in the temporary directory.
  for f in "$TMPDIR"/*.txt; do
    # Exclude header rows starting with "Grouping" and empty lines.
    # Prepend the source tree filename to each data row and append to the final file.
    grep -v "^Grouping" "$f" | grep -v "^$" | \
    awk -v t="$tree" 'BEGIN{OFS="\t"} {print t, $0}' >> "$FINAL"
  done

  # Remove the temporary directory and its contents before the next iteration.
  rm -rf "$TMPDIR"

done
