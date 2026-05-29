#!/usr/bin/env python3

# This script filters a set of genes based on their RCV (Relative Composition Variability) scores, retaining only those below the 75th percentile threshold.
# Genes with high RCV values show greater compositional heterogeneity across taxa, which can introduce systematic bias in phylogenetic inference. 
# Removing them is a common pre-filtering step to improve matrix quality.

# DEPENDENCIES:
#   numpy (pip install numpy)

# INPUT:
# A two-column tab-separated file (no header):
#   1. gene name;
#   2. RCV score (float)

# OUTPUT:
# A single-column plain-text list of gene (called list_genes_to_keep.txt) names
# whose RCV score is below the 75th percentile, one gene per line.


import numpy as np

# Input and output file paths (change with the correct names).
input_file  = "rcv_table.tsv"
output_file = "list_genes_to_keep.txt"

# Stores gene names as strings.
genes      = []
# Store corresponding RCV scores as floats.
rcv_values = []

# Read the input TSV file, skipping any empty lines.
with open(input_file) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        gene, rcv = line.split("\t")
        genes.append(gene)
        rcv_values.append(float(rcv))

# Convert the list of RCV values to a 'numpy' array for efficient computation.
rcv_array = np.array(rcv_values)

# Compute the 75th percentile of all RCV scores.
# This value is used as the upper cutoff: genes above it are considered compositionally heterogeneous and will be excluded.
threshold = np.percentile(rcv_array, 75)

# Retain only genes whose RCV score is strictly below the threshold (bottom 75%).
genes_to_keep = [g for g, r in zip(genes, rcv_array) if r < threshold]

# Write the filtered gene list to the output file, one gene name per line.
# This file can be used downstream to subset FASTA files or partition tables.
with open(output_file, "w") as out:
    for g in genes_to_keep:
        out.write(g + "\n")

print(f"Saved {len(genes_to_keep)} genes to {output_file}")
