#!/bin/bash

# Input variables
bio_h_cif=$1                    # Path to the input compressed CIF file
arpeggio_selections=$2          # Path to the arpeggio selections file
pdb_id=$3                       # PDB ID used to name output files

# Decompress the CIF file
gzip -dkc "${bio_h_cif}" > "${pdb_id}_bio-h.cif"  # Decompress the input CIF file

# Run the pdbe-arpeggio command with timeout
timeout -k 10 1s pdbe-arpeggio -sf "${arpeggio_selections}" "${pdb_id}_bio-h.cif" -i 6
exit_code=$?  # Capture the exit code of the timeout command

# Check the exit code and write to log file
if [ $exit_code -eq 124 ]; then
    # If the exit code is 124, it indicates a timeout
    echo "\"${pdb_id}\": \"timeout\"," >> "${pdb_id}_bio-h.json"
elif [ $exit_code -ne 0 ]; then
    # If the exit code is any non-zero value other than 124, it indicates a failure in pdbe-arpeggio
    echo "\"${pdb_id}\": \"arpeggio_failure\"," >> "${pdb_id}_bio-h.json"
else
    # If the exit code is 0, it indicates a successful run of pdbe-arpeggio
    echo "\"${pdb_id}\": $(cat ${pdb_id}_bio-h.json)," > "${pdb_id}_bio-h.json"

fi
#when all of the json files get combined together, need to start the combination with a semicolon, and end with a semicolon
#and potentially remove the final semicolon from the last file to be combined ? 
# Remove the decompressed CIF file
rm "${pdb_id}_bio-h.cif"  # Remove the decompressed CIF file