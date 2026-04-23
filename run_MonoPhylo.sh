# !/bin/bash

OUTDIR="/home/STUDENTI/alice.cascella/01_BUSCO_RUNs/MS80/supermatrix/monophyly_check/monophylo_results"
CHECKDIR="/home/STUDENTI/alice.cascella/01_BUSCO_RUNs/MS80/supermatrix/monophyly_check"
FINAL="$OUTDIR/final.txt"
TMPDIR="$OUTDIR/tmp_run"

# Scrivi header del file finale (una volta sola)
echo -e "source_tree\tGrouping\tNumber_Contained_Taxa\tMonophyletic\tCategory\tSupport\tNumber_Interfering_Species\tInterfering_Species" > "$FINAL"

for tree in *.treefile; do

  mkdir -p "$TMPDIR"

  python3 MonoPhylo.py \
    --tree "$CHECKDIR/$tree" \
    --out_dir "$TMPDIR" \
    --map "$CHECKDIR/Species_taxonomy_contracted3.tsv" \
    --support

  # Per ogni file di output del run, prendi solo le righe dati (non header)
  for f in "$TMPDIR"/*.txt; do
    grep -v "^Grouping" "$f" | grep -v "^$" | \
    awk -v t="$tree" 'BEGIN{OFS="\t"} {print t, $0}' >> "$FINAL"
  done

  rm -rf "$TMPDIR"

done
