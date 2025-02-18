#!/usr/bin/bash
# convert_ncbi_to_ko.sh
# This script converts a list of NCBI gene accessions to KEGG KO identifiers.

# Input file: one NCBI accession per line (without the "ncbi-geneid:" prefix)
INPUT="top_210_genes_ncbi_accessions.txt"
# Output file: tab-separated values: NCBI_accession <tab> KO_ID
OUTPUT="ncbi_to_ko.tsv"

# Write a header to the output file.
echo -e "NCBI_accession\tKO_ID" > "$OUTPUT"

# Loop over each accession in the input file.
while IFS= read -r accession; do
    # Remove any extra whitespace.
    accession=$(echo "$accession" | tr -d '[:space:]')
    
    # Skip empty lines.
    if [ -z "$accession" ]; then
        continue
    fi
    
    # Query the KEGG REST API to convert the NCBI gene accession to a KO ID.
    # The URL uses the prefix "ncbi-geneid:" for your accession.
    result=$(curl -s "https://rest.kegg.jp/link/ko/bfre:${accession}")
    
    # If a mapping is found, the output will be of the form:
    # ncbi-geneid:138633451    ko:K12345
    if [ -n "$result" ]; then
        # Extract the KO identifier (the second column) using awk.
        ko_id=$(echo "$result" | awk '{print $2}')
    else
        # If no mapping is returned, set the KO identifier as NA.
        ko_id="NA"
    fi
    
    # Write the result (NCBI accession and KO identifier) to the output file.
    echo -e "${accession}\t${ko_id}" >> "$OUTPUT"
    
    # Optional: print progress to the terminal.
    echo "Processed accession $accession -> $ko_id"
    
    # Sleep for 1 second to be courteous to the KEGG servers (adjust as needed).
    sleep 1
done < "$INPUT"

echo "Conversion complete. See the file '$OUTPUT' for results."
