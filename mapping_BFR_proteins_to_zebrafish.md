## Mapping BFR annotated proteins to protein symbols from Danio rerio

### Download Danio rerio RefSeq proteins

```{bash}
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_protein.faa.gz
```

### Make BLAST database for D. rerio proteins
```{bash}
module load ncbi-cli

crun.ncbi-cli makeblastdb -in GCF_000002035.6_GRCz11_protein.faa.gz -dbtype prot
```

Run bash script to BLAST BFR proteins against D. rerio database and extract gene symbol for top hit
```{bash}
#!/usr/bin/bash
# run_BFR_to_DRE_mapping.sh
# A fully bash-based workflow to map query gene IDs to Danio rerio gene symbols.

# ----- LOAD MODULE -----
module load ncbi-cli

# ----- CONFIGURATION -----
# Input FASTA file containing your query protein sequences.
QUERY_FASTA="GCF_029030835.1_ASM2903083v1_protein.faa"

# Local BLAST database name (assumes you ran makeblastdb on GCF_000002035.6_GRCz11_protein.faa).
BLAST_DB="./danio_rerio/GCF_000002035.6_GRCz11_protein.faa"

# Output file for BLAST results (tab-delimited format).
BLAST_OUT="BFR_to_DRE_blast_results.tsv"

# Your email address (required by NCBI).
NCBI_EMAIL="jtoy@odu.edu"

# Output file for the final mapping.
MAPPING_OUT="BFR_to_DRE_mapping.tsv"

# ----- STEP 1: Run BLASTp -----
# Run BLASTp, outputting tabular format with these columns:
# qseqid (query ID), sseqid (subject/accession), and evalue.
# We limit to the top hit with -max_target_seqs 1.
echo "Running BLASTp..."
crun.ncbi-cli blastp -query "$QUERY_FASTA" -db "$BLAST_DB" -out "$BLAST_OUT" \
       -outfmt "6 qseqid sseqid evalue" -max_target_seqs 1 -num_threads 22

# ----- STEP 2: Loop through BLAST results and fetch gene symbol -----
# Create header for mapping output.
echo -e "QueryID\tSubjectID\tEvalue\tGeneSymbol" > "$MAPPING_OUT"

echo "Processing BLAST results and mapping gene symbols..."
while IFS=$'\t' read -r qseqid sseqid evalue; do
    echo "Processing query: $qseqid, top hit: $sseqid"

    # Fetch the GenBank record for the top hit using NCBI efetch via curl.
    GENBANK=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${sseqid}&rettype=gb&retmode=text&email=${NCBI_EMAIL}")

    # Parse the GenBank record for the gene symbol.
    # We look for the first occurrence of /gene="SOMETHING".
    GENE=$(echo "$GENBANK" | grep -m1 '/gene=' | sed -n 's/.*\/gene="\([^"]*\)".*/\1/p')

    # If no gene symbol was found, set to NA.
    if [ -z "$GENE" ]; then
        GENE="NA"
    fi

    # Append the mapping to the output file.
    echo -e "${qseqid}\t${sseqid}\t${evalue}\t${GENE}" >> "$MAPPING_OUT"

    # Sleep briefly to be kind to NCBI servers (adjust if processing many queries).
    sleep 2
done < "$BLAST_OUT"

echo "Mapping complete. See ${MAPPING_OUT}."

```
