#!/usr/bin/env python3

# This script performs a monophyly check and finds the support value associated with the nodes of interest.
# It requires a tree in Newick format, on which the analysis will be performed, and a table containing every species present on the tree with the associated taxonomy. 
# 


import argparse
from ete4 import Tree

def parse_args():
    parser = argparse.ArgumentParser(
        description="Check monophyly of a clade of interest in a phylogenetic tree."
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
