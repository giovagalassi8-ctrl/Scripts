# !/usr/bin/bash

# This script automates parallel IQ-TREE3 runs across multiple supermatrices using the Dayhoff6 amino acid recoding scheme.
# Recoding reduces the 20 standard amino acids to 6 Dayhoff groups, which can mitigate compositional heterogeneity and long-branch attraction artefacts in phylogenetic inference.
# Each job is launched in a dedicated detached screen session, allowing all runs to execute simultaneously in the background.


# Set the number of cores (change based on the number of threads available).
TOTAL_THREADS=54
# Set the number of jobs (change accordingly). 
# In the example, six submatrices (folders; from 01_ to 06_) were considered into the supermatrix.
NUM_JOBS=6
THREADS_PER_JOB=$((TOTAL_THREADS / NUM_JOBS))

echo "$THREADS_PER_JOB thread used for each job"

# Set the sub-folders to process (change with your folders name).                                                                                            
FOLDERS=(01_g95 02_g90 03_g85 04_g80 05_kpic 06_kpi)

for DIR in "${FOLDERS[@]}"; do
    # Extraxt the suffix after the underscore (in this case: g95, g90, g85,...).
    SUFFIX="${DIR#*_}"
    # Set files names (MS80 refers to the supermatrix; allgenes refers to the fact that no one filters was applied to this matrix; 6aa refers to the use of the Dayhoff6 amino acidic recoding).
    # Change names accordingly.
    MATRIX="concatenated_MS80_${SUFFIX}.out"
    PARTITION="partitions_MS80_${SUFFIX}.txt"
    $OUTDIR="02_ML_MS80_${SUFFIX}_allgenes_6aa"
    PREFIX="ML_MS80_${SUFFIX}_allgenes_6aa"

    # Check if the output directory has alreay a *.treefile.
     if compgen -G "$dir/$OUTDIR/*.treefile" > /dev/null; then
        echo "Skip $DIR: there is already a treefile in $OUTDIR"
        continue
    fi

    # Creates the output directory if it does not exist.
    mkdir -p "$DIR/$OUTDIR"
    # Copy concatenated and partitions files into the output directory.
    cp "$dir/$MATRIX" "$DIR/$OUTDIR/"
    cp "$dir/$PARTITION" "$DIR/$OUTDIR/"

   # Run the iqtree command.
   CMD="iqtree3 -s $OUTDIR/$MATRIX \
        -p "$OUTDIR/$PARTITION" \
        -m Dayhoff \
        -B 1000 \
        -T $THREADS_PER_JOB \
        --prefix $PREFIX"

    echo
    echo "Run IQ-TREE for $DIR:"
    echo "$CMD"
    echo

    # Run the previous command in background, into a screen.
    screen -dmS "IQTREE_$SUFFIX" bash -c "$CMD > $OUTDIR/iqtree.log 2>&1"
done

echo "Every job has run into different screen session."
echo "Use 'screen -ls' to see every active screen session."
