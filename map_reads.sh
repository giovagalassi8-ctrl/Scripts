#!/bin/bash

# This script maps both short paired-end reads (Illumina) and long reads (PacBio/Nanopore) to a reference assembly/genome. It performs the following steps for both data types:
# 1. Map reads using minimap2 to generate a SAM file.
# 2. Convert SAM to BAM using samtools view.
# 3. Sort the BAM file.
# 4. Index the sorted BAM file.
# 5. Remove intermediate files to save space.
#
# USAGE:
# [bash] ./map_reads.sh <reference.fasta> <short_read_1.fastq> <short_read_2.fastq> <long_reads.fastq> [threads]
#
# ARGUMENTS:
# $1 : Reference genome/assembly (FASTA)
# $2 : Short Read Forward (R1)
# $3 : Short Read Reverse (R2)
# $4 : Long Reads (FASTQ or FASTA, can be .gz)
# $5 : Number of threads (Optional, default: 6)

# --- 1. Argument Check & Variable Assignment ---

if [ "$#" -lt 4 ]; then
    echo "Error: Missing arguments."
    echo "Usage: $0 <ref.fasta> <sr_1.fastq> <sr_2.fastq> <lr.fastq> [threads]"
    exit 1
fi

REF="$1"          # Reference FASTA
SR1="$2"          # Short Read R1
SR2="$3"          # Short Read R2
LR="$4"           # Long Read file
THREADS="${5:-6}" # Use 5th argument as threads, or default to 6 if not provided

# Extract the base name of the reference (e.g., "genome.fasta" -> "genome")
# This is used to name the output files consistently.
BASENAME=$(basename "$REF" | cut -d. -f1)

echo "--------------------------------------------------"
echo "Starting Mapping Pipeline"
echo "Reference: $REF"
echo "Output Prefix: $BASENAME"
echo "Using Threads: $THREADS"
echo "--------------------------------------------------"

# --- 2. Short Reads Processing ---
echo "[Short Reads] Mapping starts..."

# Map short reads using minimap2
# -ax sr : Preset for genomic short reads (SR)
# --MD   : Output the MD tag (required for some variant callers)
# -t     : Number of threads
minimap2 -ax sr --MD -t "$THREADS" "$REF" "$SR1" "$SR2" > "${BASENAME}_sr.sam"

# Convert SAM to BAM
# -S : Input is SAM (auto-detected in newer versions, but kept for compatibility)
# -b : Output is BAM
samtools view -Sb "${BASENAME}_sr.sam" > "${BASENAME}_sr.bam"

# Remove the intermediate SAM file to save disk space
rm "${BASENAME}_sr.sam"

# Sort the BAM file
# -o : Output filename
echo "[Short Reads] Sorting BAM..."
samtools sort -@ "$THREADS" -o "${BASENAME}_sr_sorted.bam" "${BASENAME}_sr.bam"

# Index the SORTED BAM file (Corrected from original script)
samtools index "${BASENAME}_sr_sorted.bam"

# Remove the unsorted BAM file
rm "${BASENAME}_sr.bam"

echo "[Short Reads] Done. Output: ${BASENAME}_sr_sorted.bam"


# --- 3. Long Reads Processing ---
echo "--------------------------------------------------"
echo "[Long Reads] Mapping starts..."

# Map long reads using minimap2
# -ax map-pb : Preset for PacBio genomic reads (use 'map-ont' for Nanopore)
minimap2 -ax map-pb --MD -t "$THREADS" "$REF" "$LR" > "${BASENAME}_lr.sam"

# Convert SAM to BAM
samtools view -Sb "${BASENAME}_lr.sam" > "${BASENAME}_lr.bam"

# Remove intermediate SAM
rm "${BASENAME}_lr.sam"

# Sort the BAM file
echo "[Long Reads] Sorting BAM..."
samtools sort -@ "$THREADS" -o "${BASENAME}_lr_sorted.bam" "${BASENAME}_lr.bam"

# Index the SORTED BAM file
samtools index "${BASENAME}_lr_sorted.bam"

# Remove the unsorted BAM file
rm "${BASENAME}_lr.bam"

echo "[Long Reads] Done. Output: ${BASENAME}_lr_sorted.bam"
echo "--------------------------------------------------"
echo "Pipeline finished successfully."
