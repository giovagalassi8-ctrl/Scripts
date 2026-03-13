#!/usr/bin/env python3

import sys
from ete4 import Tree


TREE_FILE = sys.argv[1]
TAXONOMY_FILE = sys.argv[2]
CLADE_OF_INTEREST = sys.argv[3]

t = Tree(open(TREE_FILE).read())

taxonomy = {}
with open(TAXONOMY_FILE, 'r') as f:
    next(f)  # Salta l'intestazione
    for line in f:
        parts = line.strip().split('\t')
        # Il file ha 2 colonne: specie e clade
        if len(parts) == 2 and parts[1]:  # Ignora righe senza clade
            sp, clade = parts
            taxonomy[sp] = clade

# Create a list of leaf nodes that belong to the target clade.
# Iterate over all leaves and check if their name in 'taxonomy' matches the target clade.
target_leaves = [leaf for leaf in t.leaves() if taxonomy.get(leaf.name) == CLADE_OF_INTEREST]


if len(target_leaves) >= 2:
        # Find the Most Recent Common Ancestor fot the identified leaves.
        mrca = t.common_ancestor(target_leaves)
        # Verify monophyly (check if all leaves descending from the mrca belong to the target clade).
        monophyly = all(taxonomy.get(leaf.name) == CLADE_OF_INTEREST for leaf in mrca.leaves())
        # Find the bootstrap support value from the MRCA node properties.
        bootstrap = getattr(mrca, "support", "NA")
        #Print the results
        print(f"{TREE_FILE}\t{CLADE_OF_INTEREST}\t{monophyly}\t{bootstrap}")
else:
        print(f"{TREE_FILE}\t{CLADE_OF_INTEREST}\tNot_enough_leaves\tNA")
