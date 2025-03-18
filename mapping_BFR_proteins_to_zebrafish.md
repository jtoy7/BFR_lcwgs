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


This didn't work because the .faa file excludes the genes I need to blast (those that start with LOC). New approach below:


### Download Danio rerio RefSeq proteins

```{bash}
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_cds_from_genomic.fna.gz
```

### Make BLAST database for D. rerio proteins
```{bash}
module load ncbi-cli

crun.ncbi-cli makeblastdb -in GCF_000002035.6_GRCz11_cds_from_genomic.fna -dbtype nucl
```


```
#!/usr/bin/bash
# run_BFR_to_DRE_mapping_v2.sh
# A fully bash-based workflow to map query gene IDs to Danio rerio gene symbols.


# ----- LOAD MODULE -----
module load ncbi-cli

# ----- CONFIGURATION -----
# Input file containing gene list
GENE_LIST="top_210_gene_list.txt"

# Output file
OUTPUT_FILE="BFR_to_DRE_mapping.tsv"

# Reference database for Danio rerio
BLAST_DB="./danio_rerio/GCF_000002035.6_GRCz11_cds_from_genomic.fna"

# Prepare output file
echo -e "Query_ID\tDANRE_NCBI_ID\tGene_Symbol" > $OUTPUT_FILE

# ----- Process genes -----
grep "^LOC" "$GENE_LIST" | while read -r LOC_ID; do
    echo "Processing $LOC_ID..."

    # Remove "LOC" prefix
    CLEAN_ID=${LOC_ID#LOC}

    # Fetch the nucleotide sequence from NCBI Datasets
    crun.ncbi-cli datasets download gene gene-id "$CLEAN_ID" --include gene --filename "$LOC_ID.zip"
    
    # Unzip and extract FASTA file
    unzip -o "$LOC_ID.zip" -d "$LOC_ID"
    FASTA_FILE=$(find "$LOC_ID" -name "*.fna" | head -n 1)

    if [ ! -f "$FASTA_FILE" ]; then
        echo "Warning: No sequence found for $LOC_ID"
        continue
    fi

    # Run BLAST against Danio rerio genome, keeping only the top hit
    BLAST_RESULT=$(crun.ncbi-cli blastn -query "$FASTA_FILE" -db "$BLAST_DB" -outfmt "6 qseqid sseqid ssciname scomname pident length evalue bitscore stitle" -max_target_seqs 1 -num_threads 4 | head -n 1)

    if [ -z "$BLAST_RESULT" ]; then
        echo "Warning: No BLAST hit for $LOC_ID"
        echo -e "$LOC_ID\tNo_Hit\tNo_Hit" >> $OUTPUT_FILE
        continue
    fi

    # Extract relevant fields
    NCBI_ID=$(echo "$BLAST_RESULT" | awk '{print $2}')

    # Fetch gene symbol for the BLAST hit
    GENE_SYMBOL=$(echo "$BLAST_RESULT" | sed -n 's/.*gene=\([^]]*\).*/\1/p')

    # Append to output file
    echo -e "$LOC_ID\t$NCBI_ID\t$GENE_SYMBOL" >> $OUTPUT_FILE

    # Clean up temporary files
    rm -rf "$LOC_ID" "$LOC_ID.zip"
done

echo "BLAST search completed. Results saved in $OUTPUT_FILE"
```
