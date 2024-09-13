# Load necessary libraries
library(edgeR)
library(openxlsx)
library(org.Mm.eg.db)

# Load your filtered and normalized gene data
counts <- read.csv("C:/Users/agsko/dev/pcm/filtered_genes_with_counts.csv", row.names = 1)

# Convert row names (ENSEMBL IDs) to a column
counts$ENSEMBL_ID <- rownames(counts)

# Move the ENSEMBL_ID column to the beginning
counts <- counts[, c("ENSEMBL_ID", names(counts)[-length(names(counts))])]

# Assuming 'metadata' is already loaded and contains the correct Sample ID and Condition
# The metadata should look something like this:

metadata <- data.frame(
  Sample.ID = c("X88", "X89", "X90", "X92", "X93", "X94", 
                "X97", "X99", "X101", "X103", "X104", 
                "X105", "X108", "X109", "X110"),
  Condition = c("Ure", "Ure", "Ure", 
                "Ure_DMXAA", "Ure_DMXAA", "Ure_DMXAA", 
                "Ure_RT", "Ure_RT", "Ure_RT", 
                "Ure_RT_DMXAA", "Ure_RT_DMXAA", "Ure_RT_DMXAA", 
                "PBS_RT_DMXAA", "PBS_RT_DMXAA", "PBS_RT_DMXAA")
)

# Set the row names of the metadata to match sample IDs
rownames(metadata) <- metadata$Sample.ID

# Check that the order of metadata rows matches the order of counts columns
all(rownames(metadata) == colnames(counts)[3:17])


# Create the design matrix
design <- model.matrix(~ 0 + metadata$Condition)
colnames(design) <- levels(metadata$Condition)

# View the design matrix to confirm
print(design)

# Explicitly set the column names for the design matrix
colnames(design) <- c("PBS_RT_DMXAA", "Ure", "Ure_DMXAA", "Ure_RT", "Ure_RT_DMXAA")

# Print the design matrix again to verify the column names
print(design)

# Estimate dispersion
dge <- estimateDisp(dge, design)

# You can also check the dispersion estimates
dge$common.dispersion  # Common dispersion
head(dge$tagwise.dispersion)  # Tagwise dispersion for individual genes

################Pre-processing for GOAT###################

library(org.Mm.eg.db)

# Assuming your Ensembl IDs are in the 'ENSEMBL_ID' column of dge$genes
ensembl_ids <- dge$genes$ID

# Map Ensembl IDs to Entrez IDs
entrez_ids <- mapIds(org.Mm.eg.db, keys = ensembl_ids, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")

# Map Ensembl IDs to Gene Symbols
gene_symbols <- mapIds(org.Mm.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Add these annotations to your DGEList genes data
dge$genes$EntrezID <- entrez_ids
dge$genes$GeneSymbol <- gene_symbols

# Check the first few entries to ensure the mapping was successful
head(dge$genes)

# Count the number of NA values in each column of dge$genes
na_counts <- sapply(dge$genes, function(x) sum(is.na(x)))

# Print the results
print(na_counts)

###> print(na_counts)
#ID  Gene.name   EntrezID GeneSymbol 
#0         85       9882       9882 

# Filter the rows where EntrezID is NA
missing_entrez_ids <- dge$genes[is.na(dge$genes$EntrezID), ]

# Count the number of missing Entrez IDs
num_missing_entrez <- nrow(missing_entrez_ids)

# Print the count
cat("Number of Ensembl IDs without corresponding Entrez IDs:", num_missing_entrez, "\n")

# Optionally, view the first few rows of these missing entries
head(missing_entrez_ids)

# Total number of genes
total_genes <- nrow(dge$genes)

# Number of genes without Entrez IDs
num_missing_entrez <- sum(is.na(dge$genes$EntrezID))

# Calculate the percentage of missing Entrez IDs
percentage_missing <- (num_missing_entrez / total_genes) * 100

# Print the results
cat("Total number of genes:", total_genes, "\n")
cat("Number of genes without Entrez IDs:", num_missing_entrez, "\n")
cat("Percentage of genes without Entrez IDs:", percentage_missing, "%\n")

#########> cat("Total number of genes:", total_genes, "\n")
#Total number of genes: 33268 
#> cat("Number of genes without Entrez IDs:", num_missing_entrez, "\n")
#Number of genes without Entrez IDs: 9882 
#> cat("Percentage of genes without Entrez IDs:", percentage_missing, "%\n")
#Percentage of genes without Entrez IDs: 29.70422 %

# Filter out genes without Entrez IDs
dge_filtered <- dge[!is.na(dge$genes$EntrezID), ]

# Check how many genes remain after filtering
cat("Number of genes remaining after filtering:", nrow(dge_filtered$genes), "\n")
#Number of genes remaining after filtering: 23386 

# Identify the indices of the genes that need to be updated
index_50518 <- which(dge_filtered$genes$EntrezID == "50518")
index_20997 <- which(dge_filtered$genes$EntrezID == "20997")

# Update the Gene.name and GeneSymbol for the specified Entrez IDs
dge_filtered$genes$Gene.name[index_50518] <- "Nonagouti"
dge_filtered$genes$GeneSymbol[index_50518] <- "Nonagouti"

dge_filtered$genes$Gene.name[index_20997] <- "Gm20997"
dge_filtered$genes$GeneSymbol[index_20997] <- "Gm20997"

# Check the updated entries to confirm the changes
updated_genes <- dge_filtered$genes[c(index_50518, index_20997), ]
print(updated_genes)


# Fit the GLM model
fit <- glmFit(dge, design)

# Specify the contrast you want to test (e.g., Ure_DMXAA vs Ure)
lrt <- glmLRT(fit, contrast = makeContrasts(Ure_DMXAA - Ure, levels = design))

# Extract the differential expression results
dge_results <- topTags(lrt, n = Inf)$table

# Extract the required values
effectsize <- dge_results$logFC  # log2 fold change
pvalue <- dge_results$PValue  # Unadjusted p-value
signif <- p.adjust(pvalue, method = "BH") < 0.05  # Adjusted p-value, using Benjamini-Hochberg

# Convert 'signif' to 'true'/'false' for the web tool
signif <- ifelse(signif, "true", "false")

# Ensure that the rownames of dge_results match the rownames of dge_filtered$genes
matching_indices <- match(rownames(dge_results), rownames(dge_filtered$genes))

# If necessary, reorder the rows of dge_filtered$genes to match dge_results
dge_filtered$genes <- dge_filtered$genes[matching_indices, ]

# Check the first few rownames of dge_results
head(rownames(dge_results))

# Check the first few rownames of dge_filtered$genes
head(rownames(dge_filtered$genes))

# Check if there are any genes in dge_results that are not in dge_filtered$genes
non_matching_genes <- setdiff(rownames(dge_results), rownames(dge_filtered$genes))
cat("Number of non-matching genes:", length(non_matching_genes), "\n")

# Optionally, inspect some of these non-matching genes
head(non_matching_genes)


# Combine these into a data frame with gene annotations
result_data <- data.frame(
  gene = dge$genes$EntrezID,  # Entrez ID
  symbol = dge$genes$GeneSymbol,  # Gene symbol
  effectsize = effectsize,  # log2 fold change
  pvalue = pvalue,  # Unadjusted p-value
  signif = signif  # Significance (TRUE if adjusted p-value < 0.05)
)

# Ensure there are no missing values
result_data <- na.omit(result_data)

# Check the final data frame
str(result_data)

# Check for any NA values in the data frame
missing_values_summary <- sapply(result_data, function(x) sum(is.na(x)))

# Print the summary of missing values
print(missing_values_summary)

# Save as CSV
write.csv(result_data, "Ure_DMXAAvsUre.csv", row.names = FALSE)

# List of all contrasts to analyze
contrasts_list <- list(
  Ure_DMXAAvsUre = makeContrasts(Ure_DMXAA - Ure, levels = design),
  Ure_RTvsUre = makeContrasts(Ure_RT - Ure, levels = design),
  Ure_RT_DMXAAvsUre_RT = makeContrasts(Ure_RT_DMXAA - Ure_RT, levels = design),
  Ure_RT_DMXAAvsUre = makeContrasts(Ure_RT_DMXAA - Ure, levels = design),
  UrevsPBS_RT_DMXAA = makeContrasts(Ure - PBS_RT_DMXAA, levels = design),
  Ure_DMXAAvsPBS_RT_DMXAA = makeContrasts(Ure_DMXAA - PBS_RT_DMXAA, levels = design),
  Ure_RTvsPBS_RT_DMXAA = makeContrasts(Ure_RT - PBS_RT_DMXAA, levels = design),
  Ure_RT_DMXAAvsPBS_RT_DMXAA = makeContrasts(Ure_RT_DMXAA - PBS_RT_DMXAA, levels = design),
  Ure_DMXAAvsUre_RT = makeContrasts(Ure_DMXAA - Ure_RT, levels = design),
  Ure_RT_DMXAAvsUre_DMXAA = makeContrasts(Ure_RT_DMXAA - Ure_DMXAA, levels = design)
)

# Loop over each contrast and perform the analysis
for (contrast_name in names(contrasts_list)) {
  # Get the contrast matrix
  contrast <- contrasts_list[[contrast_name]]
  
  # Perform GLM fitting and likelihood ratio test
  lrt <- glmLRT(fit, contrast = contrast)
  
  # Extract the differential expression results
  dge_results <- topTags(lrt, n = Inf)$table
  
  # Extract the required values
  effectsize <- dge_results$logFC  # log2 fold change
  pvalue <- dge_results$PValue  # Unadjusted p-value
  signif <- p.adjust(pvalue, method = "BH") < 0.05  # Adjusted p-value, using Benjamini-Hochberg
  
  # Convert 'signif' to 'true'/'false' for the web tool
  signif <- ifelse(signif, "true", "false")
  
  # Ensure that the rownames of dge_results match the rownames of dge_filtered$genes
  matching_indices <- match(rownames(dge_results), rownames(dge_filtered$genes))
  
  # If necessary, reorder the rows of dge_filtered$genes to match dge_results
  dge_filtered$genes <- dge_filtered$genes[matching_indices, ]
  
  # Combine these into a data frame with gene annotations
  result_data <- data.frame(
    gene = dge_filtered$genes$EntrezID,  # Entrez ID
    symbol = dge_filtered$genes$GeneSymbol,  # Gene symbol
    effectsize = effectsize,  # log2 fold change
    pvalue = pvalue,  # Unadjusted p-value
    signif = signif  # Significance (TRUE if adjusted p-value < 0.05)
  )
  
  # Ensure there are no missing values
  result_data <- na.omit(result_data)
  
  # Check the final data frame
  str(result_data)
  
  # Check for any NA values in the data frame
  missing_values_summary <- sapply(result_data, function(x) sum(is.na(x)))
  
  # Print the summary of missing values
  print(missing_values_summary)
  
  # Save the results as CSV
  output_filename <- paste0(contrast_name, ".csv")
  write.csv(result_data, output_filename, row.names = FALSE)
}


