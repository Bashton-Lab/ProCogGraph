#!/bin/bash

process_row() {
    bio_h_cif=$1
    arpeggio_selections=$2
    pdb_id=$3

    # decompress structure
    gzip -dkc "${bio_h_cif}" > "${pdb_id}_bio-h.cif"

    timeout -k 10 6h pdbe-arpeggio -sf "${arpeggio_selections}" "${pdb_id}_bio-h.cif" -i 6
    exit_code=$?

    #handle various exit conditions
    if [ $exit_code -eq 124 ]; then
        echo "\"${pdb_id}\": \"timeout\"" > "${pdb_id}_bio-h.json"
    elif [ $exit_code -ne 0 ]; then
        echo "\"${pdb_id}\": \"arpeggio_failure\"" > "${pdb_id}_bio-h.json"
    else
        echo "\"${pdb_id}\": $(cat ${pdb_id}_bio-h.json)" > "${pdb_id}_bio-h.json"
    fi

    # Remove the decompressed CIF file
    rm "${pdb_id}_bio-h.cif"
}

#parse the script input
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 input.csv"
    exit 1
fi

csv_file=$1

#process each row of the csv file
{
    read #skip the first row (header)
    while IFS=, read -r pdb_id assembly_id updated_mmcif bio_h_cif sifts_xml arpeggio_selections bound_entity_info; do
        process_row "$bio_h_cif" "$arpeggio_selections" "$pdb_id"
    done
} < "$csv_file"
