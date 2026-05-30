# !/bin/sh

# This script runs IQ-TREE3 with a partition model and automatic model selection (TESTNEWMERGE) on concatenated supermatrices, iterating over a list of subfolders. 
# For each subfolder, the script creates a dedicated output directory, moves the input files into it, and launches IQ-TREE3 with 1000 ultrafast bootstrap replicates.
# It requires a plain-text list of subfolder names (one per line; called subfolders_name.txt), a concatenated file and a partition file in each subfolder.
# As output returns three files:
#   *.treefile                 : best Maximum Likelihood tree;
#   *.iqtree                   : full IQ-TREE run report;
#   *.contree                  : consensus tree with bootstrap support.


# Set the matrix name (change accordingly).
matrix="MS80"

# Read subfolder names one by one from the list.
while read DIR; do
    # Enter the subfolder for this matrix/trimming combination.
    cd "$DIR"
    # Create a dedicated output directory for this IQ-TREE run.
     # (In this example: ML refers to a Maximum Likelyhood tree; allgenes refers to the fact that no one filters was applied to this matrix; PM refers to the use of Partitions Models).
    mkdir -p "ML_${matrix}_${DIR}_allgenes_PM"
    # Move the concatenated and partition files into the output directory.
    mv concatenated* "ML_${matrix}_${DIR}_allgenes_PM/"
    mv partition*    "ML_${matrix}_${DIR}_allgenes_PM/"
    # Enter the output directory before running IQ-TREE3.
    cd "ML_${matrix}_${DIR}_allspecies_PM"

    # Run IQ-TREE3:
    # -s  : input concatenated file (change accordingly).
    # -p  : input partition file (change accordingly).
    # -m TESTNEWMERGE : tests and merges partitions to find the best-fit model, reducing over-parameterisation.
    iqtree3 -s concatenated* -p partition* -m TESTNEWMERGE -T 32 -B 1000 \
        --prefix "ML_${matrix}_${DIR}_allgenes_PM"

    # Return two levels up (back to the main working directory) for the next iteration.
    cd ../../

done < subfolders_name.txt
