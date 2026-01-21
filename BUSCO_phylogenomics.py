#!/usr/bin/env python

# BUSCO_phylogenomics
# Utility script to construct species phylogenies using BUSCO results.
# Assumes the same BUSCO dataset has been used on each genome

# 2023 Jamie McGowan <jamie.mcgowan@earlham.ac.uk>
# https://github.com/jamiemcg/BUSCO_phylogenomics


import argparse
import multiprocessing as mp
from os import listdir, chdir, mkdir, system
from os.path import abspath, basename, isdir, join
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from time import gmtime, strftime

def main():
    parser = argparse.ArgumentParser(description="Perform phylogenomic reconstruction using BUSCO sequences")

    parser.add_argument("-i", "--input", type=str, help="Input directory containing completed BUSCO runs", required=True)
    parser.add_argument("-o", "--output", type=str, help="Output directory to store results", required=True)
    parser.add_argument("-t", "--threads", type=int, help="Number of threads to use", required=True)
    parser.add_argument("--supermatrix_only", help="Don't generate gene trees", action="store_true")
    parser.add_argument("--gene_trees_only", help="Don't perform supermatrix analysis", action="store_true")
    parser.add_argument("--nt", help="Align nucleotide sequences instead of amino acid sequences", action="store_true")  
    parser.add_argument("-psc", "--percent_single_copy", type=float, action="store", dest="psc", default=100.0,
                        help="BUSCO presence cut-off. BUSCOs that are complete and single-copy in at least [-psc] percent of species will be included in the contatenated alignment [default=100.0]")
    parser.add_argument("--trimal_strategy", type=str, action="store", dest="trimal_strategy", default="automated1",
                        help="trimal trimming strategy (automated1, gappyout, strict, strictplus) [default=automated1]")
    parser.add_argument("--missing_character", type=str, action="store", dest="missing_character", help="Character to represent missing data [default='?']", default="?")
    parser.add_argument("--gene_tree_program", type=str, action="store", dest="gene_tree_program", default="fasttree", help="Program to use to generate gene trees (fasttree or iqtree) [default=fasttree]>
    parser.add_argument("--busco_version_3", action="store_true", help="Flag to indicate that BUSCO version 3 was used (which has slighly different output structure)")
    
    args = parser.parse_args()

    print_message("Starting BUSCO Phylogenomics Pipeline")
    print_message("User provided arguments:", sys.argv)
    print_message("Parsed arguments:", vars(args))

    input_directory = abspath(args.input)
    working_directory = abspath(args.output)
    threads = args.threads

    # Check if the input directory exists
    if not isdir(input_directory):
        print_message("ERROR. Input BUSCO directory", input_directory, "not found")
        sys.exit()

    # Check if the output directory already exists
    if isdir(working_directory):
        print_message("ERROR. Output directory", working_directory, "already exists")
        sys.exit()
