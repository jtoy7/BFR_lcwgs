### GOWINDA setup for NCBI annotation
### Jason A. Toy
### 2025-02-10

# Load packages
library(readr)
library(dplyr)

# Set working directory
setwd("~/bfr/ncbi_annotation/")

# Read the GO annotation file (assuming tab-separated format)
go_annotations <- read_tsv("GCF_029030835.1-RS_2024_11_gene_ontology.gaf.gz", comment = "!", col_names = FALSE)

# Select relevant columns: GeneID (Column 2) and GO_ID (Column 5)
go_data <- go_annotations %>%
  dplyr::select(GeneID = X2, GO_ID = X5) %>%
  distinct()  # Remove duplicates

# Group by GO_ID and concatenate genes into space-separated format
go_data_reformatted <- go_data %>%
  group_by(GO_ID) %>%
  summarize(Gene_List = paste(unique(GeneID), collapse = " ")) %>%
  ungroup()

# Use list of GO IDs to get GO term descriptions
library(GO.db)

  # Create a helper function to get the GO term description from GO.db
  get_go_term_description <- function(go_id) {
    term_obj <- GOTERM[[go_id]]
    if (!is.null(term_obj)) {
      Term(term_obj)
    } else {
      NA_character_
    }
  }

  # Use mutate() to add description column to go_data_reformatted
  gowinda_genesetfile <- go_data_reformatted %>% 
    mutate(Description = sapply(GO_ID, get_go_term_description)) %>% # use sapply because function isn't vectorized
    dplyr::select(GO_ID, Description, Gene_List)
  

# Save reformatted associations to file
write_tsv(gowinda_genesetfile, "gowinda_genesetfile.tsv", col_names = FALSE)


# Create mapping file for original contig names to NCBI contig names from gff file
contig_name_mapping <- read_tsv("sequence_report_clean.tsv") %>% 
  dplyr::select(RefSeq_seq_accession, Sequence_name)


# Reformat GTF file to have consistent gene name format

gtf <- read_tsv("GCF_029030835.1_ASM2903083v1_genomic.gtf", comment = "#", col_names = FALSE)

# Create function to split key-value pairs in attributes column (X9) into a named list
parse_attributes_gtf <- function(attr_string) {
  attr_list <- str_split(attr_string, ";")[[1]]  # Split by semicolon
  attr_list <- str_trim(attr_list) # trims leading and trailing whitespace from a string
  attr_list <- attr_list[attr_list != ""] # remove any empty strings (which can occur from trailing semicolons).
  attr_pairs <- str_split_fixed(attr_list, "\\s+", 2)   # use "\\s+" (a regex for one or more whitespace characters) and limit to 2 pieces (key and value)
  set_names(attr_pairs[, 2], attr_pairs[, 1])  # Create named vector (key = column name, value = content)
}

# Apply function to each row and convert to a data frame
gtf_parsed <- gtf %>%
  mutate(parsed = map(X9, parse_attributes_gtf)) %>%  # Parse each attributes column into a named list stored in a new column called "parsed". Must use map() because our custom function is not vectorized (expects only a single string, not a vector of strings)
  unnest_wider(parsed) # Expand "parsed" list column into separate columns

gtf_gowinda <- gtf_parsed %>% 
  dplyr::select(1:8, db_xref, transcript_id) %>% 
  mutate(
    gene_id_clean = str_remove_all(db_xref, 'GeneID:'), # Remove the "GeneID:" prefix and any quotes from the db_xref column
    combined = paste0('gene_id ', gene_id_clean, '; transcript_id ', transcript_id, ';')     # Combine the cleaned values into the desired format
  ) %>% 
  dplyr::select(1:8, combined) %>% 
  left_join(contig_name_mapping, by = c("X1" = "RefSeq_seq_accession")) %>% # map original contig names to RefSeq names
  dplyr::select(Sequence_name, 2:9)

# Save to file
write_tsv(gtf_gowinda, "gowinda_ncbi.gtf", col_names = FALSE, escape = "none")
