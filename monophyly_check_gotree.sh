#!/bin/bash

# This script automates the monophyly check of multiple phylogenetic trees against various taxonomic groups using 'gotree'. 
# It generates a table where each row represents a specific tree, each column a taxonomic group, and the cells contains the corresponding monophyly results (true or false).

# This script must be run in a folder containing all the tree files you want to measure the monophyly.
# It also needs a folder containing .txt files, each of which refers to a particular taxonomic group present in trees. These files must contain the name of the species analyzed (with space ad a separator, instead of underscore).

# USAGE:
# [bash] run_monophyly_gotree > <output_name>


# --- INPUT FILES ---
# Stores all files ending in .treefile in the directory
trees="*.treefile"
# Stores all text files in a 'gropus' folder. In this case the folder is located one level up: change the path if necessary.
groups="../groups/*.txt"

# --- HEADER GENERATION ---

# Print the first column header ("Tree")
echo -ne "Tree"

# Loop throught each group file to create the rest of the column headers.
for g in $groups
do
    # Extract just the file name without the ".txt" extension.
    name=$(basename "$g" .txt)
    # Print a tab character followed by the group name on the same line.
    echo -ne "\t${name}"
done

# Print a clean new line to finish the header row and move to the data section.
echo ""

# --- DATA PROCESSING ---

# Loop throught each file. Each tree will represent a single row in the final table.
for t in $trees
do
    # Get the clean tree file name.
    treename=$(basename "$t")
    # Print the tree name at the beginning of the row, on the same line.
    echo -ne "${treename}"

    # The following loop test the current tree against every group file.
    for g in $groups
    do
        # Execute gotree for the monophyly check.
        result=$(gotree stats monophyletic -i "$t" -l "$g" | tail -n 1 | awk '{print $2}')
        echo -ne "\t${result}"
    done

    # After testing all groups for this specific tree, print a new line: the next tree on the loop will start on the row below.
    echo ""
done
