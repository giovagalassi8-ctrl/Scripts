# !/bin/bash

# This script trims a set of aligned FASTA files using six different ClipKIT trimming strategies in parallel, 
# writing the results into six dedicated subfolders.
# This allows downstream phylogenetic analyses to be run on multiple trimming stringencies simultaneously,
# making it easy to compare their effect on tree topology and support values.
# It requires an aligned FASTA files, one per gene, produced by an aligner such as MAFFT ( sequences/aln_*.faa).
#   The six trimming strategies applied are:
#     g95 : removes columns with >95% gaps (gappy mode);
#     g90 : removes columns with >90% gaps (gappy mode);
#     g85 : removes columns with >85% gaps (gappy mode);
#     g80 : removes columns with >80% gaps (gappy mode);
#     kpic : keeps parsimony-informative and constant sites;
#     kpi : keeps parsimony-informative sites only.
# It creates one folder per trimming strategy, each containing all trimmed alignments prefixed with "trim_".


# Path to the folder containing the aligned input FASTA files (change if needed).
INPUT_DIR="sequences"

# Write the list of output folder names to a plain-text file.
# This file can be used directly as input for downstream iteration scripts (e.g. passed to iqtree_allgenes.sh via subfolders_name.txt).
cat << EOF > folders_name.txt
01_g95
02_g90
03_g85
04_g80
05_kpic
06_kpi
EOF

# Create all six output directories.
mkdir -p 01_g95 02_g90 03_g85 04_g80 05_kpic 06_kpi

# Iterate over every aligned FASTA file in the input directory.
for FILE in ${INPUT_DIR}/aln_*.faa; do

    # Extract the filename without the directory path (e.g. aln_gene1.faa).
    BASENAME=$(basename "$FILE")

    # Apply all six trimming strategies to the current file, writing each result to the corresponding output folder with a "trim_" prefix.
    # -m : trimming mode;
    # -g : gap threshold (only used in gappy mode);
    clipkit "$FILE" -m gappy -g 0.95 -o 01_g95/trim_"$BASENAME"
    clipkit "$FILE" -m gappy -g 0.9  -o 02_g90/trim_"$BASENAME"
    clipkit "$FILE" -m gappy -g 0.85 -o 03_g85/trim_"$BASENAME"
    clipkit "$FILE" -m gappy -g 0.8  -o 04_g80/trim_"$BASENAME"
    clipkit "$FILE" -m kpic         -o 05_kpic/trim_"$BASENAME"
    clipkit "$FILE" -m kpi          -o 06_kpi/trim_"$BASENAME"

done
