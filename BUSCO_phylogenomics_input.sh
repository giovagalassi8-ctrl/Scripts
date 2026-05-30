# !/bin/bash

# Collects and organises BUSCO output files from multiple species into a single unified directory.
# After running BUSCO independently on several assemblies, each producing its own "*_busco" output folder,
# this script gathers all results into a standardised folder structure, making downstream analyses easier to manage.
# For each "*_busco" folder found in the current directory the script:
#     1. Checks whether the expected BUSCO output subfolder exists
#        (run_metazoa_odb10/busco_sequences/); if not, prints a warning.
#     2. Extracts the species name by removing the "_busco" suffix.
#     3. Creates a dedicated subdirectory for that species inside
#        all_busco_sequences/.
#     4. Copies the entire run_metazoa_odb10/ folder into it, preserving
#        the original BUSCO directory structure.


# Store the current working directory path for reference.
WORKDIR=$(pwd)

# Create the main output directory that will contain all BUSCO sequences, organised by species (change as desired).
mkdir -p all_busco_sequences

# Iterate over every folder ending in "_busco" in the current directory.
# Each folder is expected to contain the output of a single BUSCO run.
for dir in *_busco; do
    # Check whether the expected BUSCO output subfolder exists for this run.
    # If the BUSCO run was successful, the busco_sequences directory should always be present inside run_metazoa_odb10/ (change accordingly).
    if [ -d "$dir/run_metazoa_odb10/busco_sequences" ]; then

        # Extract the species name by removing the "_busco" suffix from the folder name (e.g. "Verpa_penis_busco" → "Verpa_penis").
        species=$(echo $dir | sed 's/_busco$//')
        # Create a dedicated subdirectory for this species inside the output folder.
        mkdir -p "all_busco_sequences/$species"
        # Copy the entire BUSCO results folder (run_metazoa_odb10/) into the species subdirectory, preserving the original directory structure.
        cp -r "$dir/run_metazoa_odb10" "all_busco_sequences/$species/"

    else
        # Warn the user if the expected output folder is missing, which may
        # indicate a failed or incomplete BUSCO run for this species.
        echo "WARNING: No busco_sequences in $dir"
    fi

done
