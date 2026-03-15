#!/usr/bin/env python3

# This script performs a monophyly check and finds the  bootstrap support value associated with the nodes of interest.
# It can be used considering different taxonomic ranks (those that are not present on a species tree, but which one would like to analyze).

# REQUIRED FILES:
# 1. A phylogenetic tree in Newick format;
# 2. A CSV file containing different taxonomic ranks associated to the same species that are considered on the tree (e.g. Phylum, Class, Order, Family, Genus).

# OPTIONS:
# -t, --tree : the phylogenetic tree in analysis.
# -x, --taxonomy : the CSV file containing the taxonomy associated to every species into the tree.
# -r, --rank : name of the taxonomic rank you want to evaluate (it has to correspond to the exact name of the column of interest).
# -c, --clade [optional] : name of the exactly clade you want to analyse (if not specified, the script will be run on the entire rank/column specified before).

# USAGE:
# ./monophyly_check_by_rank.py -t <TREE_NEWICK> -x <TAXONOMY_FILE> -r <COLUMN_OF_INTEREST> -c <CLADE>


# Import the argparse module to handle command-line interface arguments.
import argparse
# Import the Tree class from the ete4 library to read and analyze phylogenetic trees.
from ete4 import Tree

# Define the arguments required for the script.
def parse_args():
    parser = argparse.ArgumentParser(
        description="Check monophyly and bootstrap support of a clade of interest in a phylogenetic tree."
    )
    parser.add_argument(
        "-t", "--tree",
        required=True,
        help="Path to the tree file (Newick format)"
    )
    parser.add_argument(
        "-x", "--taxonomy",
        required=True,
        help="Path to the taxonomy file"
    )
    parser.add_argument(
        "-r", "--rank",
        required=True,
        help="Name of the taxonomic rank (header of the column) of interest"
    )
    parser.add_argument(
        "-c", "--clade",
        required=None,
        help="Name of the clade of interest"
    )
    # Parse the arguments provided by the user and return them.
    return parser.parse_args()

# Execute the parse_args function and store the result in 'args'.
args = parse_args()

TREE_FILE = args.tree
TAXONOMY_FILE = args.taxonomy

# Open the tree file, read its content, and parse it into an ete4 Tree object.
t = Tree(open(TREE_FILE).read())

# Initialize an empty dictionary to map each species to its corresponding clade.
taxonomy = {}
# Open the taxonomy file in read mode.
with open(TAXONOMY_FILE, 'r') as f:
    # Read the first line (header), remove trailing whitespace, and split by tab character.
    header = next(f).strip().split('\t')
    # Find the numeric index of the requested taxonomic rank column.
    col_index = header.index(args.rank)
    for line in f:
        # Strip trailing whitespace and split the line into parts by tab character.
        parts = line.strip().split('\t')
        # Assign the first column (index 0) as the species name.
        sp = parts[0]
        # Retrieve the clade name using the previously identified column index.
        clade = parts[col_index]
        # If the clade value is not empty, add it to the dictionary mapped to the species.
        if clade:
                taxonomy[sp] = clade

# Check if a single clade to test is specificated via the '-c' argument.
if args.clade:
    # Create a list containing only the specified clade.
    clades_to_test = [args.clade]
else:
    # Extract all unique clades from the taxonomy dictionary values to test them all.
    clades_to_test = set(taxonomy.values())

print(f"Tree_file\t{args.rank}\tMonophyletic\tBootstrap")

# Iterate over each clade that needs to be tested.
for CLADE_OF_INTEREST in clades_to_test:
    # Use a list comprehension to find all leaves in the tree that belong to the current clade.
    target_leaves = [leaf for leaf in t.leaves() if taxonomy.get(leaf.name) == CLADE_OF_INTEREST]

    # Verify if there are at least two leaves to properly define a Most Recent Common Ancestor (MRCA).
    if len(target_leaves) >= 2:
        # Find the MRCA node for the identified target leaves.
        mrca = t.common_ancestor(target_leaves)
        # Check monophyly by evaluating if all leaves under this MRCA belong exclusively to the clade of interest.
        monophyly = all(taxonomy.get(leaf.name) == CLADE_OF_INTEREST for leaf in mrca.leaves())
        # Try to retrieve the bootstrap support value of the MRCA, defaulting to "NA" if not found.
        bootstrap = getattr(mrca, "support", "NA")
        # Print the results.
        print(f"{TREE_FILE}\t{CLADE_OF_INTEREST}\t{monophyly}\t{bootstrap}")
    else:
        print(f"{TREE_FILE}\t{CLADE_OF_INTEREST}\tNot_enough_leaves\tNA")
