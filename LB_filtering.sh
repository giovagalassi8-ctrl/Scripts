# !/bin/bash

# This script calculates the LB score for a specific phylogeny, removes the outlier taxa and reconstructs the concatenated filtered alignment.
# This represent an end-to-end pipeline for long-branch (LB) detection, sequence filtering, multiple sequence alignment, trimming, and supermatrix concatenation,
# intended as a quality-control step before phylogenetic tree reconstruction.

# Run this script in a folder containig:
#    - filter_fasta_by_list_of_headers.py;
#    - LB_detect_outliers.py;
#    - A treefile on which calculate the LB score;
#    - A two-column TSV file that associates a species to a specific clade.
# Change this script accordingly to the name of your files.


# Calculates the LB scores using the long_branch_score option in Phykit. 
# Change with the correct Treefile name you want to calculate the LB score on.
phykit long_branch_score TREEFILE -v > LB_scores.tsv
# Adds the correct clade to a species based on a TSV taxonomy file (in this case called species_to_clade.tsv).
awk 'NR==FNR { clade[$1]=$2; next } 
($1 in clade) { print $1 "\t" $2 "\t" clade[$1] }
' species_to_clade.tsv LB_scores.tsv > LB_by_clade.tsv

# Run the LB_detect_outliers.py script, which detects outliers, and copy the original sequences in this folder.
python run_LB_detect_outliers.py
# Change with the correct path of the folder containing every sequence.
cp ../../sequences/*.faa .

# Filter outlier sequences and move filtered files into a new folder.
for a in *.faa; do
    python filter_fasta_by_list_of_headers.py ${a} LB_taxa_to_remove.txt > filtered_${a}
done

# Creates the alignments folder.
mkdir alignments
# Move every filtered fasta file just created into the alignments folder.
mv filtered_*.faa alignments/
# Remove every non-filtered fasta file from this folder.
rm *.faa

cd alignments
# Edit the name of every file in the alignment folder by removing the filtered_ prefix.
for filename in ./*; do mv "./$filename" "./$(echo "$filename" | sed -e 's/filtered_//g')";  done

# Align sequences with MAFFT (you can change with another aligner).
for a in *.faa; do mafft --maxiterate 1000 --localpair --thread 10 ${a} > aln_${a}; done

# Trim the just aligned sequences using clipkit. Select and uncomment the desired trimming model and comment out the rest.
# You can also add another trimming model if it is not present below.
# for a in aln_*; do clipkit  ${a} -m gappy -g 0.95; done
# for a in aln_*; do clipkit  ${a} -m gappy -g 0.90; done
# for a in aln_*; do clipkit  ${a} -m gappy -g 0.85; done
# for a in aln_*; do clipkit  ${a} -m gappy -g 0.80; done
# for a in aln_*; do clipkit  ${a} -m kpic; done
# for a in aln_*; do clipkit  ${a} -m kpi; done
# for a in aln_*; do trimal -in ${a} -out trim_${a} -gappyout; done

# Concatenate every trimmed sequences using AMAS.py .
# Change the following options:
#   -i: select the imput file based on whether ClipKIT-trimmed files (*clipkit) or trimAl-trimmed files (trim_aln*) was used in the previous step;
#       only one pattern will actually match depending on which trimming step was run above.
#   -p: change with the name of the output partitions file.
#   -t: change with the name of the output concatenated file.
AMAS.py concat -i [*clipkit][trim_aln*] -f fasta -d aa -u fasta -y nexus -p PARTITIONS.txt -t CONCATENATED.out
# Move every concatenated and partitions file into the previous folder.
mv concatenated_*.out partitions_*.txt ../
