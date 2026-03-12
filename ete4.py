import sys
from ete4 import Tree


TREE_FILE = sys.argv[1]
TAXONOMY_FILE = sys.argv[2]
CLADE_OF_INTEREST = sys.argv[3]
t = Tree(TREE_FILE)


# 'with' statement will automatically close the file after the nested block of code. 'r' option open the file in reading mode.
with open(TAXONOMY_FILE, 'r') as f:
    taxonomy = dict(line.strip().split('\t') for line in f if line.count('\t') == 2)
        # Create a dictionary to improve the research into the treefile. The main operations on a dictionary are storing a value with some key and extracting the value given the key. Kry must be unique and immutable.
        # Build the dictionary directly by filtering and splitting valid lines to store the taxonomy mapping.
        # line.strip takes the raw line, removes blank space and divide it in a list of strings using tabs as separator. This operation generates the exactly Specie-Clade pair.
# Create a list of leaf nodes that belong to the target clade.
# Iterate over all leaves and check if their name in 'taxonomy' matches the target clade.
target_leaves = [leaf for leaf in t.leaves() if taxonomy.get(leaf.name) == CLADE_OF_INTEREST]


if len(target_leaves) > 1:
        # Find the Most Recent Common Ancestor fot the identified leaves.
        mrca = t.get_common_ancestor(target_leaves)
        # Verify monophyly (check if all leaves descending from the mrca belong to the target clade).
        monophyly = all(taxonomy.get(leaf_name) == CLADE_OF_INTEREST for leaf in mrca.leaves())
        # Find the bootstrap support value from the MRCA node properties.
        bootstrap = mrca.props.get("support", "NA")
        #Print the results
        print(f"{TREE_FILE}\t{CLADE_OF_INTEREST}\t{monophyly}\t{bootstrap}")
else:
        print(f"{TREE_FILE}\t{CLADE_OF_INTEREST}\tNot_enough_leaves\tNA")
