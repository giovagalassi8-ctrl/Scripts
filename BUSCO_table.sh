# !/bin/bash

# This script parses a concatenated file of BUSCO short summary lines and converts it into a structured TSV table,
# with one row per species and one column per BUSCO category.
# This makes the results easier to inspect, filter, and use as input for downstream plotting scripts (e.g. BUSCO_plot.R).
# It requires a concatenated BUSCO short summary lines, one per species, in the format: <species> C:X%[S:X%,D:X%],F:X%,M:X%,n:X. 
# This file can be obtained by the concatenation of all the short_summary*.txt files produced after a BUSCO run for every species.

# The output file obtained is a tab-separated table with columns: Species, C, S, D, F, M, n.
# If a line does not match the expected format, all values are set to NA.


# Set the input file and the name of the output file (change accordingly).
INPUT="00_all_short_summaries.txt"
OUTPUT="01_all_short_summaries_table.tsv"

# Write the header line to the output file.
echo -e "Species\tC\tS\tD\tF\tM\tn" > "$OUTPUT"

# Parse each line of the input file with gawk and append results to the output.
awk '
{
    # The first whitespace-separated field is the species name.
    species = $1
    
    # Regex to capture all six BUSCO values from the summary string.
    # Capture groups:
    #   a[1] = C  (Complete, %)
    #   a[2] = S  (Single-copy, %)
    #   a[3] = D  (Duplicated, %)
    #   a[4] = F  (Fragmented, %)
    #   a[5] = M  (Missing, %)
    #   a[6] = n  (total BUSCOs in the dataset, integer)
    regex = "C:([0-9.]+)%\\[S:([0-9.]+)%,D:([0-9.]+)%\\],F:([0-9.]+)%,M:([0-9.]+)%,n:([0-9]+)"

    if (match($0, regex, a)) {
        # Line matched: print species name followed by the six captured values.
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", species, a[1], a[2], a[3], a[4], a[5], a[6]
    } else {
        # Line did not match the expected format: fill all value columns with NA
        # to preserve row count and flag malformed or missing entries.
        printf "%s\tNA\tNA\tNA\tNA\tNA\tNA\n", species
    }
}
' "$INPUT" >> "$OUTPUT"

echo "Table created: $OUTPUT"
