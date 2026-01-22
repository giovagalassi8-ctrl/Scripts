#!/bin/bash

output="Likelyhoood.txt"

# Svuota il file di output se esiste
> "$output"

# Cerca tutti i Base_results.txt
find . -type f -name "Base_results.txt" | while read -r file; do

    # Percorso relativo senza ./ iniziale
    rel_path="${file#./}"

    # Etichetta (00_1L/1K/2N)
    label=$(dirname "$rel_path")

    # Estrae SOLO la prima riga del file
    line=$(head -n 1 "$file")

    # Scrive nel file finale
    if [[ -n "$line" ]]; then
        echo -e "${label}\t${line}" >> "$output"
    fi

done
