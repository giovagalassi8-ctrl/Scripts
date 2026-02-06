#!/bin/bash

# ==============================================================================
# SCRIPT: check_monophyly.sh
# DESCRIZIONE:
#   Itera attraverso la struttura di cartelle (01_*/01_*, ecc.), trova il file
#   .treefile (escludendo quello _annotated) e lancia phykit monophyly_check
#   usando il file di gruppo specificato dall'utente.
#
# USO:
#   ./check_monophyly.sh <nome_file_gruppo.txt>
#   Esempio: ./check_monophyly.sh Mollusca.txt
# ==============================================================================

# 1. Controllo Argomenti
if [ "$#" -ne 1 ]; then
    echo "ERRORE: Devi specificare il nome del file di gruppo."
    echo "Uso: $0 <NomeFile.txt>"
    echo "Esempio: $0 Mollusca.txt"
    exit 1
fi

GROUP_FILENAME="$1"      # Es. Mollusca.txt
OUTPUT_FILE="report_${GROUP_FILENAME%.*}_monophyly.txt" # Es. report_Mollusca_monophyly.txt

# Intestazione del file di output
echo -e "Cartella\tRisultato_Phykit" > "$OUTPUT_FILE"

echo "--- Inizio analisi per il gruppo: $GROUP_FILENAME ---"

# 2. Iterazione attraverso le cartelle principali (01_... a 07_...)
for main_dir in 0[1-7]_*; do
    
    # Controllo che la cartella esista (per evitare errori se non ce ne sono alcune)
    [ -d "$main_dir" ] || continue

    # Iterazione attraverso le sottocartelle (01_... a 09_...) all'interno
    for sub_dir in "$main_dir"/0[0-9]_*; do
        
        [ -d "$sub_dir" ] || continue

        # 3. Trova il file .treefile corretto
        # Cerchiamo un file che finisce per .treefile ma NON per _annotated.treefile
        # head -n 1 serve per sicurezza nel caso ne trovasse più di uno, prende il primo.
        tree_file=$(find "$sub_dir" -maxdepth 1 -type f -name "*.treefile" ! -name "*_annotated.treefile" | head -n 1)

        # 4. Definisce il percorso del file dei gruppi
        group_file_path="$sub_dir/groups/$GROUP_FILENAME"

        # 5. Esegue il controllo solo se entrambi i file esistono
        if [[ -f "$tree_file" && -f "$group_file_path" ]]; then
            
            # Esegue phykit
            # Nota: phykit stampa il risultato a schermo, lo catturiamo nella variabile
            result=$(phykit monophyly_check "$tree_file" "$group_file_path")
            
            # Salva nel file di report: NomeCartella [TAB] Risultato
            echo -e "${sub_dir}\t${result}" >> "$OUTPUT_FILE"
            
            # Feedback a video (opzionale, per vedere che lavora)
            echo "[OK] Analizzato: $sub_dir"
            
        else
            echo "[SKIP] File mancanti in: $sub_dir (Tree: $tree_file | Group: $group_file_path)"
        fi

    done
done

echo "--- Completato. Risultati salvati in: $OUTPUT_FILE ---"
