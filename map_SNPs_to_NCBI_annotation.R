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



# Merge SNPs with gene data, based on contig and position criteria and accounting for gene strand
# Simultaneously classify SNP types
snps_mapped <- snps %>%
  left_join(gff_genes, by = c("refseq_contig" = "seqid")) %>%
  mutate(
    snp_type = case_when(
      pos >= start & pos <= end ~ "inside",
      strand == "+" &
        pos >= start - 200 & pos < start ~ "upstream",
      strand == "-" &
        pos > end & pos <= end + 200 ~ "upstream",
      TRUE ~ "not_in_gene"
    ),
    gene_length = end - start
  ) %>%
  filter(snp_type != "not_in_gene")


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

## 11,802 SNPs mapped to genes or 200 bp upstream of genes (11,675 inside & 127 upstream)


# Write the mapped SNP table to a file
write_tsv(snps_mapped_parsed, "snps_mapped_to_ncbi_annot_genes.tsv")


# Retain only the highest-XtX SNP for each gene to get list of unique candidate genes
snps_mapped_parsed_uniq <- snps_mapped_parsed %>%
  arrange(Dbxref, desc(XtXst)) %>%
  distinct(Dbxref, .keep_all = TRUE)

# Count number of genes with SNPs inside vs. upstream of genes
snps_mapped_parsed_uniq %>% 
  count(snp_type) %>% 
  mutate(percent = 100 * n / sum(n))

## 5,678 SNPs mapped to genes or 200 bp upstream of genes (5,578 inside & 100 upstream)


# Count genes with at least one intragenic SNP, irrespective of which SNP had the highest XtX
snps_mapped_parsed %>%
  group_by(Dbxref) %>%
  summarise(
    has_inside = any(snp_type == "inside"),
    has_upstream = any(snp_type == "upstream"),
    .groups = "drop"
  ) %>%
  summarise(
    total_genes = n(),
    genes_with_inside = sum(has_inside),
    genes_only_upstream = sum(!has_inside & has_upstream),
    genes_only_inside = sum(has_inside & !has_upstream),
    genes_with_both = sum(has_inside & has_upstream),
    pct_with_inside = 100 * mean(has_inside)
  )

# total_genes genes_with_inside genes_only_upstream genes_only_inside genes_with_both pct_with_inside
# <int>             <int>               <int>             <int>           <int>           <dbl>
#  5678              5602                  76              5554              48            98.7

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



### A few additional investigations of candidate snp to gene mapping

# Check for SNPs mapping to more than one gene
# ...before thinning to one SNP per gene
snps_mapped_parsed %>%
  group_by(refseq_contig, pos) %>%
  summarise(
    n_genes = n_distinct(Dbxref),
    .groups = "drop"
    ) %>%
  filter(n_genes > 1)

# A tibble: 492 × 3
# refseq_contig      pos n_genes
# <chr>            <dbl>   <int>
#   1 NW_027210600.1 9614722       2
# 2 NW_027210600.1 9621431       2
# 3 NW_027210601.1  868778       2
# 4 NW_027210601.1  891223       2
# 5 NW_027210601.1  907190       2
# 6 NW_027210601.1  908672       2
# 7 NW_027210601.1  909889       2
# 8 NW_027210601.1  934197       2
# 9 NW_027210601.1 7238109       3
# 10 NW_027210601.1 7238112       3
# # ℹ 482 more rows

# Before thinning to one SNP per gene, there are 493 SNPs that map to more than one gene

# ...after thinning to one SNP per gene
snps_mapped_parsed_uniq %>%
  group_by(refseq_contig, pos) %>%
  summarise(
    n_genes = n_distinct(Dbxref),
    .groups = "drop"
    ) %>%
  filter(n_genes > 1)

# A tibble: 123 × 3
# refseq_contig       pos n_genes
# <chr>             <dbl>   <int>
# 1 NW_027210600.1  9614722       2
# 2 NW_027210601.1   908672       2
# 3 NW_027210601.1  7352090       2
# 4 NW_027210601.1  8011300       2
# 5 NW_027210601.1 11740491       2
# 6 NW_027210602.1  3343976       2
# 7 NW_027210602.1  5159906       2
# 8 NW_027210602.1  5759981       2
# 9 NW_027210603.1  1987578       2
# 10 NW_027210603.1  4568346       2
# # ℹ 113 more rows

# After thinning to one SNP per gene, there are 123 SNPs that map to more than one gene


# Check the position of these mulitply mapping SNPs in each associated gene
snps_mapped_parsed_uniq %>%
  group_by(refseq_contig, pos) %>%
  filter(n_distinct(Dbxref) > 1) %>%
  summarise(
    n_genes = n_distinct(Dbxref),
    mapping_types = paste(
      sort(unique(snp_type)),
      collapse = " + "
      ),
    .groups = "drop"
    ) %>%
  count(n_genes, mapping_types)

# A tibble: 5 × 3
# n_genes mapping_types         n
#     <int> <chr>             <int>
# 1       2 inside              107
# 2       2 inside + upstream    11
# 3       2 upstream              2
# 4       3 inside                2
# 5       3 inside + upstream     1

# After thinning to one SNP per gene, most of the 123 multiply mapping SNPs are intragenic in both genes


# Check how many of these overlap genes are protein coding
snps_mapped_parsed_uniq %>%
  group_by(refseq_contig, pos) %>%
  filter(n_distinct(Dbxref) > 1) %>%
  summarise(
    n_genes = n_distinct(Dbxref),
    biotypes = paste(sort(gene_biotype), collapse = " + "),
    .groups = "drop"
    ) %>%
  count(n_genes, biotypes, sort = TRUE)

# A tibble: 6 × 3
#     n_genes biotypes                                             n
#      <int> <chr>                                            <int>
# 1       2 protein_coding + protein_coding                     83
# 2       2 lncRNA + protein_coding                             35
# 3       3 protein_coding + protein_coding + protein_coding     2
# 4       2 protein_coding + snoRNA                              1
# 5       2 protein_coding + tRNA                                1
# 6       3 protein_coding + protein_coding + tRNA               1

# The vast majority are overlaps between two protein coding genes, and almost all the rest are overlaps between a protein coding gene and an lncRNA


# Check which multiply mapping SNPs, and how many, are represented in the list of top 50 candidate genes
snps_mapped_parsed_uniq %>%
  #filter(gene_biotype == "protein_coding") %>%
  arrange(desc(XtXst)) %>%
  slice_head(n = 51) %>% 
  group_by(refseq_contig, pos) %>%
  filter(n_distinct(Dbxref) > 1) %>%
  arrange(desc(XtXst)) %>%
  select(Name, Dbxref, refseq_contig, pos, XtXst, strand, start, end, snp_type)

# A tibble: 10 × 9
# Groups:   refseq_contig, pos [5]
#   Name         Dbxref           refseq_contig      pos XtXst strand   start     end snp_type
#   <chr>        <chr>            <chr>            <dbl> <dbl> <chr>    <dbl>   <dbl> <chr>   
# 1 tln1         GeneID:138637122 NW_027211032.1  101465  90.3 -        54753  147390 inside  
# 2 LOC138637126 GeneID:138637126 NW_027211032.1  101465  90.3 +         5610  127392 inside  
# 3 erbb2        GeneID:138616989 NW_027210610.1 4002572  76.3 +      3988176 4013651 inside  
# 4 LOC138617029 GeneID:138617029 NW_027210610.1 4002572  76.3 -      3951876 4198893 inside  
# 5 dock1        GeneID:138629424 NW_027210714.1  542736  75.2 +       461781  635134 inside  
# 6 LOC138629430 GeneID:138629430 NW_027210714.1  542736  75.2 -       536286  566191 inside  
# 7 syngap1b     GeneID:138630594 NW_027210736.1  315779  71.2 +       220744  349854 inside  
# 8 LOC138630626 GeneID:138630626 NW_027210736.1  315779  71.2 -       146755  350162 inside  
# 9 cog5         GeneID:138633807 NW_027210836.1   67639  70.7 -        10462   88808 inside  
# 10 gpr22a       GeneID:138633810 NW_027210836.1   67639  70.7 +        66575   72772 inside 

# Five SNPs are represented by two genes each in the list of top 50 genes
# In all cases but one, there is one assigned gene symbol and one unassigned gene
# In all cases, the SNP maps inside both assigned genes and genes are on opposite strands


# Check which multiply mapping SNPs, and how many, are represented in the list of top 207 candidate genes (XtX >= 60)
snps_mapped_parsed_uniq %>%
  #filter(gene_biotype == "protein_coding") %>%
  arrange(desc(XtXst)) %>%
  filter(XtXst >= 60) %>% 
  group_by(refseq_contig, pos) %>%
  filter(n_distinct(Dbxref) > 1) %>%
  arrange(desc(XtXst)) %>%
  select(Name, Dbxref, refseq_contig, pos, XtXst, strand, start, end, snp_type) %>% 
  print(n = 30)

# A tibble: 24 × 9
# Groups:   refseq_contig, pos [12]
#   Name         Dbxref           refseq_contig      pos XtXst strand   start     end snp_type
#   <chr>        <chr>            <chr>            <dbl> <dbl> <chr>    <dbl>   <dbl> <chr>   
# 1 tln1         GeneID:138637122 NW_027211032.1  101465  90.3 -        54753  147390 inside  
# 2 LOC138637126 GeneID:138637126 NW_027211032.1  101465  90.3 +         5610  127392 inside  
# 3 erbb2        GeneID:138616989 NW_027210610.1 4002572  76.3 +      3988176 4013651 inside  
# 4 LOC138617029 GeneID:138617029 NW_027210610.1 4002572  76.3 -      3951876 4198893 inside  
# 5 dock1        GeneID:138629424 NW_027210714.1  542736  75.2 +       461781  635134 inside  
# 6 LOC138629430 GeneID:138629430 NW_027210714.1  542736  75.2 -       536286  566191 inside  
# 7 syngap1b     GeneID:138630594 NW_027210736.1  315779  71.2 +       220744  349854 inside  
# 8 LOC138630626 GeneID:138630626 NW_027210736.1  315779  71.2 -       146755  350162 inside  
# 9 cog5         GeneID:138633807 NW_027210836.1   67639  70.7 -        10462   88808 inside  
# 10 gpr22a       GeneID:138633810 NW_027210836.1   67639  70.7 +        66575   72772 inside  
# 11 arel1        GeneID:138632737 NW_027210792.1  435914  68.4 +       436088  447912 upstream
# 12 fcf1         GeneID:138632739 NW_027210792.1  435914  68.4 -       434047  435958 inside  
# 13 sgo2         GeneID:138634244 NW_027210854.1   56568  64.9 +        51040   57180 inside  
# 14 LOC138634248 GeneID:138634248 NW_027210854.1   56568  64.9 +        44118   64100 inside  
# 15 upf2         GeneID:138636130 NW_027210954.1  229793  64.7 +       227425  243761 inside  
# 16 LOC138636131 GeneID:138636131 NW_027210954.1  229793  64.7 -       121023  251453 inside  
# 17 srrm4        GeneID:138631415 NW_027210755.1  412385  63.3 +       396638  656419 inside  
# 18 taok3a       GeneID:138631416 NW_027210755.1  412385  63.3 -       382436  458746 inside  
# 19 LOC138616752 GeneID:138616752 NW_027210610.1 3617604  62.9 +      3604356 3624867 inside  
# 20 vkorc1       GeneID:138616974 NW_027210610.1 3617604  62.9 +      3616385 3619559 inside  
# 21 tshz1        GeneID:138620556 NW_027210601.1 8011300  62.5 -      7982930 8013258 inside  
# 22 LOC138620573 GeneID:138620573 NW_027210601.1 8011300  62.5 +      7837525 8128533 inside  
# 23 pou6f2       GeneID:138619707 NW_027210620.1 1293057  61.5 +      1267281 1334128 inside  
# 24 LOC138619713 GeneID:138619713 NW_027210620.1 1293057  61.5 +      1228488 1387284 inside 

# 12 SNPs map to two genes each, so the top 207 genes are represented by 195 unique SNPs
# Only one of these genes has a "upstream" mapping
# Three SNPs map to two overlapping genes on the same strand. In each case, there is one gene with an assigned symbol and one unnamed gene.




### Analyze and plot relationship between SNP number and gene length

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
  
  # dispersion ratio is 1.361 which is greater than 1 and p-value < 0.001 so data is overdispersed. Poisson model not a good fit.

  
# try negative binomial model
  library(MASS)
  
  nb_model <- glm.nb(snp_count ~ gene_length, data = snp_counts_per_gene)
  summary(nb_model)
  
  # check for overdispersion
  check_overdispersion(nb_model)
  
# dispersion ratio = 1.556, p-value = 0.024; More overdispersed than Poisson model
  
# compare models
  AIC(p_model, nb_model)

  # df              AIC
  # p_model   2 19379.1372209469
  # nb_model  3 18953.8723014400
  
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
  check_overdispersion(root_model) # no overdispersion detected! (ratio = 1.001, p-value = 0.952)
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
  geom_point(aes(x = gene_length, y = snp_count), alpha = 0.3) +
  labs(title = "SNP Count vs Gene Length", x = "Gene Length", y = "SNP Count") +
  theme_bw()


# With model and labels
library(ggrepel)

ggplot(snp_counts_per_gene) +
  geom_point(aes(x = gene_length, y = snp_count, color = outlier), alpha = 0.2) +
  scale_color_manual(values = c("black", "red")) +  # Black for non-outliers, red for outliers
  geom_smooth(aes(x = gene_length, y = snp_count)) +
  geom_line(aes(x = gene_length, y = y_pred), color = "blue", linewidth = 1, linetype = "dashed") + # Model fit
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
  ylim(c(0, 0.0025))
  
  
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

