#!/bin/bash

trees="*.treefile"
groups="../groups/*.txt"

# header
echo -ne "Tree"

for g in $groups
do
    name=$(basename "$g" .txt)
    echo -ne "\t${name}"
done

echo ""

# righe
for t in $trees
do
    treename=$(basename "$t")
    echo -ne "${treename}"

    for g in $groups
    do
        result=$(gotree stats monophyletic -i "$t" -l "$g" | tail -n 1 | awk '{print $2}')
        echo -ne "\t${result}"
    done

    echo ""
done
