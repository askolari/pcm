# Load necessary libraries
library(edgeR)

# Assuming dge_filtered is your filtered DGEList object from the previous step

# Step 1: Calculate normalized counts (CPM) for the filtered genes
cpm_filtered <- cpm(dge_filtered, log = FALSE)  # log = FALSE returns raw CPMs

# Step 2: Recalculate log2CPM for the filtered genes (optional, if you need it)
log2cpm_filtered <- cpm(dge_filtered, log = TRUE, prior.count = 1)

# Step 3: Create a data frame with gene IDs, names, and normalized counts
filtered_gene_data_with_norm_counts <- data.frame(
  ID = rownames(dge_filtered),
  Gene.name = dge_filtered$genes$Gene.name,
  cpm_filtered  # This adds the filtered normalized CPM counts as columns
)

# Step 4: Save the filtered gene data with normalized counts to a CSV file
write.csv(filtered_gene_data_with_norm_counts, "filtered_genes_with_normalized_counts_test.csv", row.names = FALSE)

# Optional: Save the log2CPM values to a CSV file
filtered_gene_data_with_log2cpm <- data.frame(
  ID = rownames(dge_filtered),
  Gene.name = dge_filtered$genes$Gene.name,
  log2cpm_filtered  # This adds the log2 CPM values as columns
)
write.csv(filtered_gene_data_with_log2cpm, "filtered_genes_with_log2cpm.csv", row.names = FALSE)

# Assuming `dge_filtered` contains the differential expression results from edgeR


# Load edgeR library
library(edgeR)
library(openxlsx)

# Example: Load your count data and metadata (adjust paths as necessary)
counts <- read.csv("C:/Users/agsko/dev/pcm/filtered_genes_with_normalized_counts.csv", row.names = 1)  # Ensure the first column is the row names (gene IDs)
metadata <- read.xlsx("C:/Users/agsko/dev/pcm/coldata_WES_RNAseq.xlsx", sheet=1)  # Sample ID and Condition

# Convert row names (ENSEMBL IDs) to a column
counts$ENSEMBL_ID <- rownames(counts)

# Move the ENSEMBL_ID column to the beginning
counts <- counts[, c("ENSEMBL_ID", names(counts)[-length(names(counts))])]

# Check the structure to confirm the changes
str(counts)

# Create the DGEList object with the normalized data
dge <- DGEList(counts = counts[, 3:17], genes = counts[, 1:2])

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


# Make sure that the row names of metadata correspond to the column names of the counts
rownames(metadata) <- metadata$Sample.ID  # Adjust to match your column names
metadata <- metadata[ , -1]  # Remove the Sample ID column since it's now the row names

# Estimate dispersion
dge <- estimateDisp(dge, design)

# You can also check the dispersion estimates
dge$common.dispersion  # Common dispersion
head(dge$tagwise.dispersion)  # Tagwise dispersion for individual genes


library(org.Mm.eg.db)

# Assuming your Ensembl IDs are in the 'ENSEMBL_ID' column of dge$genes
ensembl_ids <- dge$genes$ENSEMBL_ID

# Map Ensembl IDs to Entrez IDs
entrez_ids <- mapIds(org.Mm.eg.db, keys = ensembl_ids, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")

# Map Ensembl IDs to Gene Symbols
gene_symbols <- mapIds(org.Mm.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Add these annotations to your DGEList genes data
dge$genes$EntrezID <- entrez_ids
dge$genes$GeneSymbol <- gene_symbols

# Check the first few entries to ensure the mapping was successful
head(dge$genes)

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

# Save as CSV
write.csv(result_data, "Ure_DMXAAvsUre.csv", row.names = FALSE)

####Ure vs Ure_RT

# Specify the contrast
lrt <- glmLRT(fit, contrast = makeContrasts(Ure_RT - Ure, levels = design))

# Extract the differential expression results
dge_results <- topTags(lrt, n = Inf)$table

# Extract the required values
effectsize <- dge_results$logFC
pvalue <- dge_results$PValue
signif <- p.adjust(pvalue, method = "BH") < 0.05

signif <- ifelse(signif, "true", "false")

result_data <- data.frame(
  gene = dge$genes$EntrezID,
  symbol = dge$genes$GeneSymbol,
  effectsize = effectsize,
  pvalue = pvalue,
  signif = signif
)

result_data <- na.omit(result_data)

write.csv(result_data, "Ure_RTvsUre.csv", row.names = FALSE)

######Ure vs Ure_RT_DMXAA
# Specify the contrast
lrt <- glmLRT(fit, contrast = makeContrasts(Ure_RT_DMXAA - Ure, levels = design))

# Extract the differential expression results
dge_results <- topTags(lrt, n = Inf)$table

# Extract the required values
effectsize <- dge_results$logFC
pvalue <- dge_results$PValue
signif <- p.adjust(pvalue, method = "BH") < 0.05

signif <- ifelse(signif, "true", "false")

result_data <- data.frame(
  gene = dge$genes$EntrezID,
  symbol = dge$genes$GeneSymbol,
  effectsize = effectsize,
  pvalue = pvalue,
  signif = signif
)

result_data <- na.omit(result_data)

write.csv(result_data, "Ure_RT_DMXAAvsUre.csv", row.names = FALSE)

##### Ure vs PBS_RT_DMXAA
# Specify the contrast
lrt <- glmLRT(fit, contrast = makeContrasts(PBS_RT_DMXAA - Ure, levels = design))

# Extract the differential expression results
dge_results <- topTags(lrt, n = Inf)$table

# Extract the required values
effectsize <- dge_results$logFC
pvalue <- dge_results$PValue
signif <- p.adjust(pvalue, method = "BH") < 0.05

signif <- ifelse(signif, "true", "false")

result_data <- data.frame(
  gene = dge$genes$EntrezID,
  symbol = dge$genes$GeneSymbol,
  effectsize = effectsize,
  pvalue = pvalue,
  signif = signif
)

result_data <- na.omit(result_data)

write.csv(result_data, "PBS_RT_DMXAAvsUre.csv", row.names = FALSE)

#######Ure_DMXAA vs Ure_RT
# Specify the contrast
lrt <- glmLRT(fit, contrast = makeContrasts(Ure_RT - Ure_DMXAA, levels = design))

# Extract the differential expression results
dge_results <- topTags(lrt, n = Inf)$table

# Extract the required values
effectsize <- dge_results$logFC
pvalue <- dge_results$PValue
signif <- p.adjust(pvalue, method = "BH") < 0.05

signif <- ifelse(signif, "true", "false")

result_data <- data.frame(
  gene = dge$genes$EntrezID,
  symbol = dge$genes$GeneSymbol,
  effectsize = effectsize,
  pvalue = pvalue,
  signif = signif
)

result_data <- na.omit(result_data)

write.csv(result_data, "Ure_RTvsUre_DMXAA.csv", row.names = FALSE)

###########Ure_DMXAA vs Ure_RT_DMXAA
# Specify the contrast
lrt <- glmLRT(fit, contrast = makeContrasts(Ure_RT_DMXAA - Ure_DMXAA, levels = design))

# Extract the differential expression results
dge_results <- topTags(lrt, n = Inf)$table

# Extract the required values
effectsize <- dge_results$logFC
pvalue <- dge_results$PValue
signif <- p.adjust(pvalue, method = "BH") < 0.05

signif <- ifelse(signif, "true", "false")

result_data <- data.frame(
  gene = dge$genes$EntrezID,
  symbol = dge$genes$GeneSymbol,
  effectsize = effectsize,
  pvalue = pvalue,
  signif = signif
)

result_data <- na.omit(result_data)

write.csv(result_data, "Ure_RT_DMXAAvsUre_DMXAA.csv", row.names = FALSE)

#####Ure_DMXAA vs PBS_RT_DMXAA
# Specify the contrast
lrt <- glmLRT(fit, contrast = makeContrasts(PBS_RT_DMXAA - Ure_DMXAA, levels = design))

# Extract the differential expression results
dge_results <- topTags(lrt, n = Inf)$table

# Extract the required values
effectsize <- dge_results$logFC
pvalue <- dge_results$PValue
signif <- p.adjust(pvalue, method = "BH") < 0.05

signif <- ifelse(signif, "true", "false")

result_data <- data.frame(
  gene = dge$genes$EntrezID,
  symbol = dge$genes$GeneSymbol,
  effectsize = effectsize,
  pvalue = pvalue,
  signif = signif
)

result_data <- na.omit(result_data)

write.csv(result_data, "PBS_RT_DMXAAvsUre_DMXAA.csv", row.names = FALSE)

####Ure_RT vs Ure_RT_DMXAA
# Specify the contrast
lrt <- glmLRT(fit, contrast = makeContrasts(Ure_RT_DMXAA - Ure_RT, levels = design))

# Extract the differential expression results
dge_results <- topTags(lrt, n = Inf)$table

# Extract the required values
effectsize <- dge_results$logFC
pvalue <- dge_results$PValue
signif <- p.adjust(pvalue, method = "BH") < 0.05

signif <- ifelse(signif, "true", "false")

result_data <- data.frame(
  gene = dge$genes$EntrezID,
  symbol = dge$genes$GeneSymbol,
  effectsize = effectsize,
  pvalue = pvalue,
  signif = signif
)

result_data <- na.omit(result_data)

write.csv(result_data, "Ure_RT_DMXAAvsUre_RT.csv", row.names = FALSE)

#########Ure_RT vs PBS_RT_DMXAA
# Specify the contrast
lrt <- glmLRT(fit, contrast = makeContrasts(PBS_RT_DMXAA - Ure_RT, levels = design))

# Extract the differential expression results
dge_results <- topTags(lrt, n = Inf)$table

# Extract the required values
effectsize <- dge_results$logFC
pvalue <- dge_results$PValue
signif <- p.adjust(pvalue, method = "BH") < 0.05

signif <- ifelse(signif, "true", "false")

result_data <- data.frame(
  gene = dge$genes$EntrezID,
  symbol = dge$genes$GeneSymbol,
  effectsize = effectsize,
  pvalue = pvalue,
  signif = signif
)

result_data <- na.omit(result_data)

write.csv(result_data, "PBS_RT_DMXAAvsUre_RT.csv", row.names = FALSE)

##########Ure_RT_DMXAA vs PBS_RT_DMXAA
# Specify the contrast
lrt <- glmLRT(fit, contrast = makeContrasts(PBS_RT_DMXAA - Ure_RT_DMXAA, levels = design))

# Extract the differential expression results
dge_results <- topTags(lrt, n = Inf)$table

# Extract the required values
effectsize <- dge_results$logFC
pvalue <- dge_results$PValue
signif <- p.adjust(pvalue, method = "BH") < 0.05

signif <- ifelse(signif, "true", "false")

result_data <- data.frame(
  gene = dge$genes$EntrezID,
  symbol = dge$genes$GeneSymbol,
  effectsize = effectsize,
  pvalue = pvalue,
  signif = signif
)

result_data <- na.omit(result_data)

write.csv(result_data, "PBS_RT_DMXAAvsUre_RT_DMXAA.csv", row.names = FALSE)
