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
# -c, --clade [optional] : name of the exactly clade you want to analyse (if not specified, the script will be run on the entire rank/column specified before)

# USAGE:
# ./monophyly_check_by_rank.py -t <TREE_NEWICK> -x <TAXONOMY_FILE> -r <COLUMN_OF_INTEREST> -c <CLADE>


import argparse
from ete4 import Tree

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
    return parser.parse_args()

args = parse_args()
TREE_FILE = args.tree
TAXONOMY_FILE = args.taxonomy

t = Tree(open(TREE_FILE).read())

taxonomy = {}
with open(TAXONOMY_FILE, 'r') as f:
    header = next(f).strip().split('\t')
    col_index = header.index(args.rank)
    for line in f:
        parts = line.strip().split('\t')
        sp = parts[0]
        clade = parts[col_index]
        if clade:
                taxonomy[sp] = clade

if args.clade:
    clades_to_test = [args.clade]
else:
    clades_to_test = set(taxonomy.values())

print(f"Tree_file\t{args.rank}\tMonophyletic\tBootstrap")

for CLADE_OF_INTEREST in clades_to_test:
    target_leaves = [leaf for leaf in t.leaves() if taxonomy.get(leaf.name) == CLADE_OF_INTEREST]

    if len(target_leaves) >= 2:
        mrca = t.common_ancestor(target_leaves)
        monophyly = all(taxonomy.get(leaf.name) == CLADE_OF_INTEREST for leaf in mrca.leaves())
        bootstrap = getattr(mrca, "support", "NA")
        print(f"{TREE_FILE}\t{CLADE_OF_INTEREST}\t{monophyly}\t{bootstrap}")
    else:
        print(f"{TREE_FILE}\t{CLADE_OF_INTEREST}\tNot_enough_leaves\tNA")
