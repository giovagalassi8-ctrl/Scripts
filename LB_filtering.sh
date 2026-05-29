#!/bin/bash


# Calculates the LB scores using the long_branch_score option in Phykit. Change with the correct Treefile name.
phykit long_branch_score TREEFILE -v > LB_scores.tsv

# Adds the correct clade to a species based on a TSV taxonomy file (in this case called species_to_clade.tsv).
awk 'NR==FNR { clade[$1]=$2; next } ($1 in clade) { print $1 "\t" $2 "\t" clade[$1] }' species_to_clade.tsv LB_scores.tsv > LB_by_clade.tsv

# Run the LB_detect_outliers.py script, which detects outliers, and copy the original sequences in this folder.
python run_LB_detect_outliers.py

# Change with the correct path of the folder containing every sequence.
cp ../../sequences/*.faa .

# Filter outlier sequences, move filtered files into a new folder.
for a in *.faa; do
    python run_LB_filter_fasta_by_list_of_headers.py ${a} 03_LB_taxa_to_remove.txt > filtered_${a}
done

mkdir alignments
mv filtered_*.faa alignments/
rm *.faa
cd alignments
# edit filename
for filename in ./*; do mv "./$filename" "./$(echo "$filename" | sed -e 's/filtered_//g')";  done

# allinea <<<<<<<<<<<<<<<<<<<<<<<<<< MODIFICA --thread se serve
for a in *.faa; do mafft --maxiterate 1000 --localpair --thread 10 ${a} > aln_${a}; done

# trimma secondo il modello che serve: <<<<<<<<<<<<<<<<<<<<<<<<< MODIFICA!!
#for a in aln_*; do clipkit  ${a} -m gappy -g 0.95; done
#for a in aln_*; do clipkit  ${a} -m gappy -g 0.90; done
#for a in aln_*; do clipkit  ${a} -m gappy -g 0.85; done
#for a in aln_*; do clipkit  ${a} -m gappy -g 0.80; done
#for a in aln_*; do clipkit  ${a} -m kpic; done
for a in aln_*; do clipkit  ${a} -m kpi; done
# for a in aln_*; do trimal -in ${a} -out trim_${a} -gappyout; done

# concatena: cambia il nome dei file di output <<<<<<<<<<<<<<<<<<<<<<<<<< MODIFICA -i, -p e -t!
AMAS.py concat -i [*clipkit][trim_aln*] -f fasta -d aa -u fasta -y nexus -p partitions_MS90_kpi_lb.txt -t concatenated_MS90_kpi_lb.out
mv concatenated_*.out partitions_*.txt ..

# Fine!
