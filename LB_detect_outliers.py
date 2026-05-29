#!/usr/bin/env python3

# This script identifies and remove long-branch (LB) outlier taxa from a phylogenetic datase, 
# working indipendently whithin each taxonomic clade.
# Long branches can distort phylogenetic inference, so flagging and removing them before
# tree reconstruction is a common quality-control step.

# The outlier detection strategy adapts to clade size:
#   - Clades with 4+ taxa : any taxon with an LB score above the 75th
#                           percentile of positive LB values is flagged.
#   - Clades with 2-3 taxa: only the single taxon with the highest positive
#                           LB score is flagged.
#   - Clades with 1 taxon : skipped (no comparison possible).
# In all cases, taxa with an LB score ≤ 0 are excluded from consideration.

# DEPENDENCIES:
#  numpy (pip install numpy)

# INPUT:
#  Tab-separated file (LB_by_clade.tsv) with NO header, three columns:
#    1. taxon name
#    2. LB score (float)
#    3. clade name

# OUTPUT:
#   - Tab-separated file (LB_outliers.tsv) with no header, that lists all flagged outliers: taxon, LB score, clade.
#   - A single-column plain-text (LB_taxa_to_remove.txt) list of outlier taxon names, ready to be used as an exclusion list (e.g. as input to filter_fasta_by_list_of_headers.py).


import numpy as np
from collections import defaultdict

# Input file: LB scores per taxon, grouped by clade (no header, tab-separated).
LB_BY_CLADE_FILE = "LB_by_clade.tsv"

# Output file: full outlier records (taxon, LB score, clade), no header.
OUTLIERS_FILE = "LB_outliers.tsv"

# Output file: plain list of outlier taxon names (one per line), for downstream filtering.
TAXA_REMOVE_FILE = "LB_taxa_to_remove.txt"

# The following object will hold tuples of (taxon, lb_score, clade) for every row in the input file.
data = []

with open(LB_BY_CLADE_FILE) as f:
    for line in f:
        # Strip trailing newline/spaces and split on tab into exactly three fields.
        taxon, lb, clade = line.strip().split("\t")
        # Cast lb to float; keep taxon and clade as strings.
        data.append((taxon, float(lb), clade))

# defaultdict(list) automatically creates an empty list for any new clade key, avoiding manual initialisation checks.
by_clade = defaultdict(list)

for taxon, lb, clade in data:
    by_clade[clade].append((taxon, lb))

# The following object accumulates flagged (taxon, lb, clade) tuples.
outliers = []

for clade, entries in by_clade.items():
    n_taxa = len(entries)  # Total number of taxa in this clade (including LB ≤ 0).

    # Only consider taxa with a strictly positive LB score.
    # Zero or negative values indicate no long-branch issue and are excluded.
    positive_entries = [(taxon, lb) for taxon, lb in entries if lb > 0]

    if n_taxa >= 4:
        # For clades large enough for a meaningful statistical threshold, use the 75th percentile (Q3) of positive LB scores as the cutoff.
        # Any taxon exceeding Q3 is considered an outlier.
        if positive_entries:
            pos_lbs = np.array([lb for _, lb in positive_entries])
            q3 = np.percentile(pos_lbs, 75)  # 75th percentile threshold.

            for taxon, lb in positive_entries:
                if lb > q3:
                    outliers.append((taxon, lb, clade))

    elif n_taxa in [2, 3]:
        # For very small clades, a percentile threshold is not meaningful.
        # Instead, only the single taxon with the highest positive LB score is flagged, to avoid removing too large a proportion of the clade.
        if positive_entries:
            max_taxon, max_lb = max(positive_entries, key=lambda x: x[1])
            outliers.append((max_taxon, max_lb, clade))

# Clades with only 1 taxon are implicitly skipped: no within-clade comparison is possible.


# Full records for inspection or downstream analysis; no header line.
with open(OUTLIERS_FILE, "w") as out:
    for taxon, lb, clade in outliers:
        out.write(f"{taxon}\t{lb}\t{clade}\n")

# Single-column list of taxon names to exclude.
# This file can be passed directly to filter_fasta_by_list_of_headers.py
# to remove the corresponding sequences from a FASTA alignment.
with open(TAXA_REMOVE_FILE, "w") as out:
    for taxon, _, _ in outliers:
        out.write(f"{taxon}\n")


# Write a summary in standard output.
print(f"Outlier intra-clade scritti in: {OUTLIERS_FILE}")
print(f"Lista taxa da rimuovere scritta in: {TAXA_REMOVE_FILE}")
print(f"Totale outlier individuati: {len(outliers)}")
