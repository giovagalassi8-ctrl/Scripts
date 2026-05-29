# !/usr/bin/bash

# This script automates parallel IQ-TREE3 runs across multiple trimmed supermatrices, each representing a different alignment-trimming stringency or method.
# Each job is launched in a dedicated detached screen session, allowing all runs to execute simultaneously in the background.


# Set the number of cores (change based on the number of threads available).
TOTAL_THREADS=54
# Set the number of jobs (change accordingly). 
# In the example, six submatrices (folders; from 01_ to 06_) were considered into the supermatrix.
NUM_JOBS=6
THREADS_PER_JOB=$((TOTAL_THREADS / NUM_JOBS))

# Restricts the substitution model search to LG only, skipping the full ModelFinder scan. 
# Use a comma-separated list (e.g. "LG,WAG,JTT") to allow multiple models.
MODELS="LG"

echo "$THREADS_PER_JOB thread used for each job"
echo "Allowed Models: $MODELS"
echo

# Set the sub-folders to process (change with your folders name).                                                                                            
FOLDERS=(01_g95 02_g90 03_g85 04_g80 05_kpic 06_kpi)

for DIR in "${FOLDERS[@]}"; do
    # Extraxt the suffix after the underscore (in this case: g95, g90, g85,...).
    SUFFIX="${DIR#*_}"
    
    # Set the output directory and if does not exist, it is created. Change the name as you desire. 
    # (In this example: ML refers to a Maximum Likelyhood tree; MS80 refers to a supermatrix; allgenes refers to the fact that no one filters was applied to this matrix; PM refers to the use of a Partitions Model).
    OUTDIR="$DIR/01_ML_MS80_${SUFFIX}_allgenes_PM"
    mkdir -p "$OUTDIR"

    # Check if the output directory has alreay a *.treefile.
    if compgen -G "$OUTDIR/*.treefile" > /dev/null; then
        echo "Skip $DIR: there is already a treefile in $OUTDIR"
        continue
    fi

    # Set the input file.
    # Change accordingly with the correct name of yours concatenated and partitions files. 
    MAT="$DIR/concatenated_${SUFFIX}.out"
    PART="$DIR/partitions_${SUFFIX}.txt"

    # Check the existence of these files.
    if [[ ! -f "$MAT" ]]; then
        echo "Matrice mancante: $MAT  (skip $DIR)"
        continue
    fi
    if [[ ! -f "$PART" ]]; then
        echo "File partizioni mancante: $PART  (skip $DIR)"
        continue
    fi

    # Copy these files into the output directory.
    cp "$MAT" "$OUTDIR/"
    cp "$PART" "$OUTDIR/"

    # Set the final treefile name. Change as desired (in this example, PM refers to a partitions model)
    PREFIX="$OUTDIR/ML_MS80_${SUFFIX}_PM"

    # Run the iqtree command.
    CMD="iqtree3 -s $OUTDIR/$(basename $MAT) \
        -p $OUTDIR/$(basename $PART) \
        -m MFP -mset $MODELS \
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
