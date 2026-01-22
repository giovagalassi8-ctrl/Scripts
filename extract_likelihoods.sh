#!/bin/bash

# Nome del file di output
output="Likelyhoood.txt"

# Se il file esiste già, lo svuota per evitare duplicati
> "$output"

# Cerca ricorsivamente tutti i file chiamati Base_results.txt
find . -type f -name "Base_results.txt" | while read -r file; do

    # Estrae il percorso relativo senza ./ iniziale
    rel_path="${file#./}"

    # Estrae la cartella fino a 00_xL/xK/xN
    label=$(dirname "$rel_path")

    # Prende la prima riga che contiene la stringa richiesta
    line=$(grep -m 1 "Model Base Final Likelyhood (-lnL)" "$file")

    # Se la riga esiste, la scrive nel file di output con l'etichetta
    if [[ -n "$line" ]]; then
        echo -e "${label}\t${line}" >> "$output"
    fi

done
