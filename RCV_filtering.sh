# !/bin/bash

# This script filters the dataset removing the genes with high RCV using the RCV_keep75.sh script.
# This represent a pipeline for RCV (Relative Composition Variability) based gene filtering, followed by multiple sequence alignment, trimming, and supermatrix concatenation. 
# Genes with high RCV scores show greater compositional heterogeneity across taxa, which can introduce systematic bias in phylogenetic inference; removing them improves matrix quality.
# It needs the same two-column tab-separated file (no header) containing gene name and RCV score required by RCV_keep75.py script.


# Creates a directory to contain the alignments.
mkdir alignments_rcv
# Run the RCV_keep75.py script to filters a set of genes based on their RCV.
python RCV_keep75.py
# Remove every 'trim_aln_' from the names of the genes to keep into the text file.
sed -i 's/trim_aln_//g' list_genes_to_keep.txt

# Copy every genes into the list_of_genes_to_keep.txt file into the alignment folder previously created.
# Change with the correct path of the directory containing the sequence files (in this example, they are into the sequence/ folder).
while read i; do cp ../../sequences/${i} alignments_rcv/; done < list_genes_to_keep.txt

cd alignments_rcv
# Align sequences with MAFFT (you can change with another aligner).
for a in *.faa; do mafft --maxiterate 1000 --localpair --thread 20 ${a} > aln_${a}; done

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
