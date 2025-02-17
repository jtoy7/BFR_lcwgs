rm(list = ls())

# Load libraries
library(tidyverse)

# Set working directory
setwd("~/bfr/ncbi_annotation/")

# Load SNP file and GFF file
snps <- read_tsv("xtx_fpca_common_outlier_snp_list_refseq_contigs.tsv", col_names = TRUE)
gff_genes <- read_tsv("GCF_029030835.1_ASM2903083v1_genomic_genes.gff", col_names = FALSE)

# Ensure proper column names for the GFF file
colnames(gff_genes) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")


# Create a function to determine SNP type
determine_snp_type <- function(snp_pos, gene_start, gene_end) {
  if (snp_pos >= gene_start & snp_pos <= gene_end) {
    return("inside")
  } else if (snp_pos >= gene_start - 200 & snp_pos < gene_start) {
    return("upstream")
  } else {
    return("not_in_gene")
  }
}

# Merge SNPs with gene data based on contig and position criteria
snps_mapped <- snps %>%
  left_join(gff_genes, by = c("refseq_contig" = "seqid")) %>%
  filter(pos >= start - 200 & pos <= end) %>%
  mutate(
    snp_type = mapply(determine_snp_type, pos, start, end),
    gene_length = end - start
  ) 

# Create function to split key-value pairs into a named list
parse_attributes <- function(attr_string) {
  attr_list <- str_split(attr_string, ";")[[1]]  # Split by semicolon
  attr_pairs <- str_split(attr_list, "=", n = 2, simplify = TRUE)  # Split by =
  set_names(attr_pairs[, 2], attr_pairs[, 1])  # Create named vector (key = column name, value = content)
}

# Apply function to each row and convert to a data frame
snps_mapped_parsed <- snps_mapped %>%
  mutate(parsed = map(attributes, parse_attributes)) %>%  # Parse each attributes column into a named list stored in a new column called "parsed". Must use map() because our custom function is not vectorized (expects only a single string, not a vector of strings)
  unnest_wider(parsed) %>%  # Expand "parsed" list column into separate columns
  dplyr::select(original_contig, refseq_contig, pos, snp_type, Name, description, gene_length, everything())

# Count number of SNPs inside vs. upstream of genes
snps_mapped_parsed %>% count(snp_type)

## 11,814 SNPs mapped to genes or 200 bp upstream of genes (11,675 inside & 139 upstream)


# Write the mapped SNP table to a file
write_tsv(snps_mapped_parsed, "snps_mapped_to_ncbi_annot_genes.tsv")


# Remove consecutive duplicates (multiple SNPs in the same gene annotation) to get list of unique candidate genes
snps_mapped_parsed_uniq <- snps_mapped_parsed %>%
  arrange(Dbxref, desc(XtXst)) %>%  # Ensure data frame is sorted by refseq_contig, then pos
  mutate(prev_geneid = lag(Dbxref)) %>%  # Create a new column with the previous row's value of "Name"
  filter(Dbxref != prev_geneid | is.na(prev_geneid)) %>%  # Keep only the first row of consecutive duplicates
  dplyr::select(-prev_geneid)  # Remove the helper column

# Count number of genes with SNPs inside vs. upstream of genes
snps_mapped_parsed_uniq %>% count(snp_type)

## 5,691 SNPs mapped to genes or 200 bp upstream of genes (5,575 inside & 116 upstream)


# Write the deduplicated mapped SNP table to a file
write_tsv(snps_mapped_parsed_uniq, "snps_mapped_to_ncbi_annot_genes_uniq.tsv")


# Write a version of the deduplicated mapped SNP table that is sorted by XtXst to a file
snps_mapped_parsed_uniq_XtXst_sorted <- snps_mapped_parsed_uniq %>% arrange(desc(XtXst))
write_tsv(snps_mapped_parsed_uniq_XtXst_sorted, "snps_mapped_to_ncbi_annot_genes_uniq_XtXst_sorted.tsv")

# Write a list of unique candidate genes to a file
uniq_gene_list <- snps_mapped_parsed_uniq$Name
write(uniq_gene_list, file = "uniq_candidate_genes_ncbi_annot.txt", sep = "\n")

# Write two-column table of candidate gene lengths to a file
gene_len_table <- snps_mapped_parsed_uniq %>% dplyr::select(Name, gene_length)
write_tsv(gene_len_table, "uniq_candidate_genes_ncbi_annot_with_lengths.tsv")



### Analyze and plot relationsip between SNP number and gene length

# Summarize data to count SNPs per gene
snp_counts_per_gene <- snps_mapped_parsed %>% 
  group_by(Dbxref, Name, description, gene_length) %>% 
  summarise(snp_count = n())

# Plot the distribution of the response variable (SNP count) and look at the shape
ggplot(snp_counts_per_gene) +
  geom_histogram(aes(x = snp_count), binwidth = 1)

  # Poisson or Negative binomial distribution could fit the data

# Create fit for data using a Poisson model
p_model <- glm(snp_count ~ gene_length, family = poisson, data = snp_counts_per_gene) # Create a poisson model because snp_count is a count variable
  
# check for overdispersion
  mean_snp <- mean(snp_counts_per_gene$snp_count)
  var_snp <- var(snp_counts_per_gene$snp_count)
  cat("Mean SNP count:", mean_snp, "\nVariance SNP count:", var_snp)

  # check if residual deviance is much greater than degrees of freedom 
  summary(p_model) # residual deviance = 5636.3, df = 5689; since residual deviance ≈ df, overdispersion does NOT appear severe

  # formal test for overdispersion
  library(performance)
  check_overdispersion(p_model)  # Replace 'model' with your GLM object
  
  # dispersion ratio is 1.363 which is greater than 1 and p-value < 0.001 so data is overdispersed. Poisson model not a good fit.

  
# try negative binomial model
  library(MASS)
  
  nb_model <- glm.nb(snp_count ~ gene_length, data = snp_counts_per_gene)
  summary(nb_model)
  
  # check for overdispersion
  check_overdispersion(nb_model)
  
  # dispersion ratio = 1.567, p-value = 0.008; More overdispersed than Poisson model
  
# compare models
  AIC(p_model, nb_model)

  #         df      AIC
  #p_model   2 19417.37
  #nb_model  3 18988.74
  
  # simulate and visually check residuals of both models
  library(DHARMa)
  # For Poisson
  poisson_res <- simulateResiduals(fittedModel = p_model, plot = TRUE)
  
  # For Negative Binomial
  nb_res <- simulateResiduals(fittedModel = nb_model, plot = TRUE)
  
  
  # General Additive Model (GAM)
  library(mgcv)
  gam_model <- gam(snp_count ~ s(gene_length), family = nb, data = snp_counts_per_gene)
  summary(gam_model)
  plot(gam_model, se = TRUE, shade = TRUE, main = "Smooth Relationship") # The y-axis (s(gene_length)) shows the contribution of the smooth term for gene_length to the predicted value of the response variable (snp_count).
  gam_res <- simulateResiduals(gam_model, plot = TRUE) # not working
  
  # Root function model
  root_model <- glm(snp_count ~ I(sqrt(gene_length)), data = snp_counts_per_gene)
  summary(root_model)
  check_overdispersion(root_model) # no overdispersion detected! (ratio = 1.001)
  root_res <- simulateResiduals(root_model, plot = TRUE)
  
  # linear model
  l_model <- lm(snp_count ~ gene_length, data = snp_counts_per_gene)
  summary(l_model)
  l_res <- simulateResiduals(l_model, plot = TRUE)

  
  # Compare models with AIC
  AIC(l_model, p_model, nb_model, gam_model, root_model) # gam model has lowest AIC, then NB, then poisson, root is highest
  
  # Compare model deviances
  anova(l_model, p_model, nb_model, gam_model, root_model, test = "Chisq")
  
  
    
# recreate dataframe to include model-predicted values and calculate residuals    
snp_counts_per_gene <- snps_mapped_parsed %>% 
  group_by(Dbxref, Name, description, gene_length) %>% 
  summarise(snp_count = n()) %>% 
  ungroup() %>% # ungroup data so it doesn't interfere with mutate
  mutate(
    snps_per_gene_normalized = (snp_count/gene_length),
    y_pred = predict(nb_model, type = "response"),
    residuals = residuals(nb_model, type = "pearson"),  # calculate standardized residuals
    outlier = residuals > 3  # Flag points more than 3 SD above predicted
    )

# Plot just the data
ggplot(snp_counts_per_gene) +
  geom_point(aes(x = gene_length, y = snp_count, alpha = 0.3)) +
  labs(title = "SNP Count vs Gene Length", x = "Gene Length", y = "SNP Count") +
  theme_bw()


# With model and labels

ggplot(snp_counts_per_gene) +
  geom_point(aes(x = gene_length, y = snp_count, color = outlier), alpha = 0.2) +
  scale_color_manual(values = c("black", "red")) +  # Black for non-outliers, red for outliers
  geom_line(aes(x = gene_length, y = y_pred), color = "blue", linewidth = 1) + # Model fit
  #geom_text(data = filter(snp_counts_per_gene, outlier==TRUE), aes(x = gene_length, y = snp_count, label = Name), size = 2, color = "red") +  # Labels for outliers
  geom_text_repel(
    data = filter(snp_counts_per_gene, outlier == TRUE) %>% arrange(desc(residuals)) %>% slice_head(n = 30),  # Only label first 30 outliers
    aes(x = gene_length, y = snp_count, label = if_else(startsWith(Name, "LOC"), description, Name)),  # Use description if Name starts with "LOC"
    size = 3, 
    color = "darkred", 
    #box.padding = 0.5,  # Space around text
    #point.padding = 0.2,  # Space from point
    #nudge_y = 2,  # Move text upward slightly
    max.overlaps = Inf  # Avoid dropping labels
  ) +
  ylim(c(-1,50)) +
  labs(title = "SNP Count vs Gene Length 3SD Outliers", x = "Gene Length", y = "SNP Count") +
  theme_bw()

# write 3 SD outlier snps to file
write_tsv(snp_counts_per_gene %>% filter(outlier==TRUE), "snp_count_per_gene_outliers_3sd.tsv")

# create 3 SD outlier object
threesd_outliers <- snp_counts_per_gene %>% filter(outlier==TRUE)

# count how many 3sd outlier genes are in the top 200 of the list of candidate genes (ranked by XtXst)
threesd_outliers$Dbxref %in% head(snps_mapped_parsed_uniq_XtXst_sorted$Dbxref, n =200) %>% sum()
# count = 23 3sd outlier genes in top 200 of all candidate genes

# which ones are they?
top200xtx_genes <- head(snps_mapped_parsed_uniq_XtXst_sorted$Dbxref, n = 200)

common_genes <- intersect(threesd_outliers$Dbxref, top200xtx_genes)

threesd_top200xtx_overlap <- snps_mapped_parsed_uniq_XtXst_sorted %>% 
  filter(Dbxref %in% common_genes)


# Plot SNP count per gene normalized for gene length

ggplot(snp_counts_per_gene) +
  geom_histogram(aes(x = snps_per_gene_normalized), binwidth = .000001)+
  xlim(c(0,0.0005))
  
ggplot(snp_counts_per_gene) +
  geom_point(aes(x = gene_length, y = snps_per_gene_normalized), alpha = 0.2) +
  ylim(c(0, 0.0025)) +
  
  
# Try to fit a hyperbolic model
h_model <- nls(snps_per_gene_normalized ~ a / (gene_length + b), start = list(a = 1, b = 1), data = snp_counts_per_gene)
summary(h_model)

normalized <- snp_counts_per_gene %>% 
  mutate(
    norm_pred = predict(h_model, type = "response"),
    norm_residuals = residuals(h_model, type = "pearson"),  # calculate standardized residuals
    norm_outlier = norm_residuals > 3  # Flag points more than 3 SD above predicted
  ) %>% 
  arrange(desc(norm_outlier), desc(norm_residuals))

ggplot(normalized) +
  geom_point(aes(x = gene_length, y = snps_per_gene_normalized, color = norm_outlier), alpha = 0.2) +
  scale_color_manual(values = c("black", "red")) +  # Black for non-outliers, red for outliers
  geom_line(aes(x = gene_length, y = norm_pred), color = "blue") +
  geom_text_repel(
    data = filter(normalized, norm_outlier == TRUE) %>% arrange(desc(residuals)) %>% slice_head(n = 30),  # Only label first 30 outliers
    aes(x = gene_length, y = snps_per_gene_normalized, label = if_else(startsWith(Name, "LOC"), description, Name)),  # Use description if Name starts with "LOC"
    size = 3, 
    color = "darkred", 
    #box.padding = 0.5,  # Space around text
    #point.padding = 0.2,  # Space from point
    #nudge_y = 2,  # Move text upward slightly
    max.overlaps = Inf  # Avoid dropping labels
  ) +
  #ylim(c(0, 0.0025)) +
  labs(title = "Normalized SNP Count vs Gene Length 3SD Outliers", x = "Gene Length", y = "SNP Count") +
  theme_bw()





## Functional Enrichment Analysis with GOseq

# load GO associations from NCBI annotation

library(dplyr)

gaf <- read_tsv("GCF_029030835.1-RS_2024_11_gene_ontology.gaf.gz", comment = "!", col_names = FALSE)  # ignores lines starting with "!"

gene2go <- gaf %>% 
  dplyr::select(gene_id = X2, go_term = X5) %>% 
  distinct() %>% 
  group_by(gene_id) %>% 
  summarise(go_terms = paste(go_term, collapse = ","), .groups = 'drop')


# write reformatted go associations to file
write_tsv(gene2go, "gene2go.txt", col_names = FALSE)

# convert to list format expected by GOseq
gene2go_list <- setNames(strsplit(gene2go$go_terms, split = ","), gene2go$gene_id)

# Check the structure; it should be a list.
str(gene2go_list)

# create list of gene IDs for candidate genes
candidate_ids <- snps_mapped_parsed_uniq %>% 
  dplyr::select(Dbxref) %>% 
  mutate(gene_id = str_split_i(Dbxref, ":", 2)) %>% 
  dplyr::select(gene_id)

# create list of gene IDs for all genes
all_ids <- gff_genes %>% 
  mutate(parsed = map(attributes, parse_attributes)) %>%  # Parse each attributes column into a named list stored in a new column called "parsed". Must use map() because our custom function is not vectorized (expects only a single string, not a vector of strings)
  unnest_wider(parsed) %>%  # Expand "parsed" list column into separate columns
  dplyr::select(Dbxref) %>% 
  mutate(gene_id = str_split_i(Dbxref, ":", 2)) %>% 
  dplyr::select(gene_id) %>% 
  arrange(gene_id)

# create logical vector of candidate genes
gene_vector <- as.integer(all_ids$gene_id %in% candidate_ids$gene_id)
table(gene_vector) # check counts
names(gene_vector) <- all_ids$gene_id


# created named vector of gene lengths
all_gene_lengths <- gff_genes %>% 
  mutate(parsed = map(attributes, parse_attributes)) %>%  # Parse each attributes column into a named list stored in a new column called "parsed". Must use map() because our custom function is not vectorized (expects only a single string, not a vector of strings)
  unnest_wider(parsed) %>%  # Expand "parsed" list column into separate columns
  mutate(gene_id = str_split_i(Dbxref, ":", 2),
         length = end - start) %>% 
  arrange(gene_id) %>% 
  dplyr::select(gene_id, length) %>% 
  deframe() #converts two-column df to a named vector
  




# compute the probability weighting function for gene length (estimates the probability that a gene is selected as a candidate as a function of gene length)
library(goseq)

pwf <- nullp(gene_vector, bias.data=all_gene_lengths)


# run enrichment analysis with goseq() command
goseq_results <- goseq(pwf, gene2cat=gene2go_list)

# Adjust p-values using the Benjamini-Hochberg method (FDR)
goseq_results$padj <- p.adjust(goseq_results$over_represented_pvalue, method="BH")

# write results to file
write_tsv(goseq_results, "GOseq_enrichment_analysis_results.tsv", col_names = TRUE)


# Display significant GO categories (e.g., padj < 0.05)
significant_GO <- subset(goseq_results, padj < 0.05)
  # no significantly enriched GO categories after FDR correction

