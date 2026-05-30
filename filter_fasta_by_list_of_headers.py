# !/usr/bin/env python3

# Script source: https://bioinformatics.stackexchange.com/questions/3931/remove-delete-sequences-by-id-from-multifasta

# This script filters a multi-FASTA file by removing sequences whose headers match entries in a user-provided exclusion list (that you can obtain using LB_detect_outliers.py). 
# Sequences NOT in the list are written to stdout, making it easy to redirect output to a new file.
# The output is a filtered FASTA sequences printed to stdout; redirect with ">" to save to a file.
# A warning is printed to stderr if any IDs in the exclusion list were not found in the input FASTA (e.g. due to typos or mismatched headers).

# DEPENDENCIES:
#   Biopython (pip install biopython)

# USAGE:
# python filter_fasta_by_list_of_headers.py input.fasta exclusion_list.txt > filtered.fasta

# ARGUMENTS:
#   input.fasta          - Multi-FASTA file containing all sequences to filter
#   exclusion_list.txt   - Plain-text file with one sequence ID per line (no ">" prefix)


# Biopython library for parsing biological sequence file formats.
from Bio import SeqIO
import sys

# Parse the input FASTA file (sys.argv[1]) as an iterator.
# SeqIO.parse yields SeqRecord objects one at a time, which is memory-efficient for large files since it avoids loading all sequences at once.
ffile = SeqIO.parse(sys.argv[1], "fasta")
# Build a set of sequence IDs to exclude from the output.
header_set = set(line.strip() for line in open(sys.argv[2]))

# Iterate over each sequence record in the FASTA file.
for seq_record in ffile:
    try:
        # Attempt to remove the current sequence's name from the exclusion set.
        # seq_record.name corresponds to the first word of the FASTA header (after ">").
        # If the name IS in the set, it gets removed and we skip printing (i.e. the sequence is filtered out).
        # No output is produced for excluded sequences.
        header_set.remove(seq_record.name)
    except KeyError:
        # KeyError means the sequence ID was NOT in the exclusion set,
        # so this sequence should be kept: print it in FASTA format to stdout.
        print(seq_record.format("fasta"))
        continue

# After processing all sequences, check whether any IDs from the exclusion list were never matched.
# A non-empty header_set indicates unmatched IDs, which may signal typos, version mismatches, or IDs absent from the input FASTA.
# The warning is sent to stderr so it does not contaminate the FASTA output on stdout.
if len(header_set) != 0:
    print(
        len(header_set),
        'of the headers from list were not identified in the input fasta file.',
        file=sys.stderr
    )
