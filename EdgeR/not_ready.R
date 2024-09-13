getwd() #should be pcm
setwd("C:/Users/agsko/dev/pcm/EdgeR")

# Load the SummarizedExperiment object from the .RDS file
sexp <- readRDS("C:/Users/agsko/dev/pcm/SummarizedExperiment_rnaseq_with_metadata.RDS")

# Load necessary libraries
library(edgeR)

# Extract gene metadata from the SummarizedExperiment object
all.gene.meta <- as.data.frame(rowData(sexp))  # Extract rowData (gene metadata)

# Create a DGEList object from the raw counts
dge <- DGEList(counts = assay(sexp, "counts"))

# Add gene metadata to the DGEList object
dge$genes <- all.gene.meta

#Normalisation

# Remove genes with zero counts across all samples
keep <- rowSums(dge$counts) > 0
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Calculate normalization factors to account for differences in library size
dge <- calcNormFactors(dge, method = "TMM")

# Calculate the average log2 CPM values for each gene

a <- aveLogCPM(dge)

library(ggplot2)

# Set the threshold for presence
avelogcpm_presence_threshold <- -1

# Create ECDF data
ecdf_data <- data.frame(
  x = sort(a),
  y = ecdf(a)(sort(a))
)

# Plotting both Histogram and ECDF with corrected code
p <- list(
  Histogram = ggplot(data.frame(logCPM = a)) +
    aes(x = logCPM) +
    geom_histogram(aes(y = 100 * (..count..) / sum(..count..)), binwidth = 0.25, boundary = 0) +
    geom_vline(xintercept = avelogcpm_presence_threshold, color = "red", linetype = "dashed") +
    xlab("Average log2CPM") + ylab("Percent of genes in bin") +
    coord_cartesian(xlim = quantile(a, c(0, 0.995)), ylim = c(0, 25)) +
    labs(title = "Average gene Log2CPM distribution for genes with at least 1 read") +
    theme(plot.caption = element_text(hjust = 0)),
  
  ECDF = ggplot(ecdf_data) +
    aes(x = x, y = y * 100) +
    geom_step() +
    geom_vline(xintercept = avelogcpm_presence_threshold, color = "red", linetype = "dashed") +
    xlab("Average logCPM") + ylab("Percent of genes with smaller average logCPM") +
    coord_cartesian(xlim = quantile(a, c(0, 0.995))) +
    labs(title = "Empirical Cumulative Distribution Function of gene LogCPM values",
         subtitle = "for genes with at least 1 read") +
    theme(plot.caption = element_text(hjust = 0))
)

# Print the plots
lapply(p, print)

# Step 1: Filter genes based on the threshold
present_genes <- rownames(dge)[a >= avelogcpm_presence_threshold]

# Subset the DGEList object to keep only the present genes
dge_filtered <- dge[present_genes, , keep.lib.sizes=FALSE]

################IMPORTANT##################

# Save the dge_filtered object to an RDS file
saveRDS(dge_filtered, "dge_filtered_normalized.RDS")

# Step 2: Extract the corresponding gene IDs and names
filtered_gene_data <- data.frame(
  ID = rownames(dge_filtered),
  Gene.name = dge_filtered$genes$Gene.name
)
#now extract gene IDs, names and raw counts
filtered_gene_data_with_counts <- data.frame(
  ID = rownames(dge_filtered),
  Gene.name = dge_filtered$genes$Gene.name,
  dge_filtered$counts  # This adds the filtered raw counts as columns
)

# Step 3: Save the filtered gene data to a CSV file
write.csv(filtered_gene_data, "filtered_genes_only.csv", row.names = FALSE)

write.csv(filtered_gene_data_with_counts, "filtered_genes_with_counts.csv", row.names = FALSE)

#Data exploration based on multidimensional scaling
# Assuming 'metadata' contains the sample conditions
# Create a color vector based on the condition in your metadata
sample_conditions <- as.factor(metadata$Condition)

# Create a named vector mapping conditions to colors
condition_colors <- c(
  "Ure" = "blue",
  "Ure_DMXAA" = "red",
  "Ure_RT" = "green",
  "Ure_RT_DMXAA" = "purple",
  "PBS_RT_DMXAA" = "orange"
)

# Use the match function to ensure correct color mapping
colors <- condition_colors[match(sample_conditions, names(condition_colors))]

# MDS plot for DGE after TMM normalization (all counts) with colors for conditions
plotMDS(dge, main = "MDS Plot After TMM Normalisation", cex = 0.8, col = colors)

# MDS plot for DGE after filtering and normalization with colors for conditions
plotMDS(dge_filtered, main = "MDS Plot After Filtering", cex = 0.8, col = colors)

# Optionally add a legend to the plot to indicate which color corresponds to which condition
legend("topright", legend = unique(sample_conditions), col = unique(colors), pch = 16, title = "Condition")


# MDS plot for DGE after TMM normalization (all counts)
plotMDS(dge, cex=0.7, main = "MDS Plot After TMM Normalisation")

# MDS plot for DGE after filtering and normalization
plotMDS(dge_filtered, main = "MDS Plot After Filtering")

# Check variance in expression for all samples
expression_variance <- apply(dge_filtered$counts, 2, var)

# Print variance of X94 compared to X92 and X93
print(expression_variance[c("X92", "X93", "X94")])

#X92      X93      X94 
#26394813 24203102 18167424 

# Perform PCA on normalized counts (after filtering)
pca_result <- prcomp(t(dge_filtered$counts), scale. = TRUE)

# Create a data frame for PCA results
pca_data <- as.data.frame(pca_result$x)
pca_data$Sample <- rownames(pca_data)

# Add the sample conditions to the PCA data
pca_data$Condition <- metadata$Condition

percent_var <- round(100 * summary(pca_result)$importance[2, 1:2], 2) # for PC1 and PC2

# Plot PCA and color by conditions
ggplot(pca_data, aes(x = PC1, y = PC2, label = Sample, color = Condition)) +
  geom_point(aes(color=Condition)) +
  geom_text(vjust = 0.5, hjust = -0.2, show.legend = FALSE) +
  theme_minimal() +
  ggtitle("PCA Plot of Samples Colored by Condition") +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +  # Add variance percentage to x-axis
  ylab(paste0("PC2: ", percent_var[2], "% variance"))  

# Extract the eigenvalues (sdev^2 gives the variance explained)
pc_eigenvalues <- data.frame(PC = 1:length(pca_result$sdev),
                             eigenvalue = pca_result$sdev^2)
# Calculate the percentage variance explained by each PC
pc_eigenvalues <- pc_eigenvalues %>%
  mutate(pct = eigenvalue / sum(eigenvalue) * 100)  # Fraction of variance explained

# Calculate cumulative percentage variance explained
pc_eigenvalues <- pc_eigenvalues %>%
  mutate(pct_cum = cumsum(pct))  # Cumulative variance explained

# Plot the variance explained by each PC
library(ggplot2)
ggplot(pc_eigenvalues, aes(x = PC)) +
  geom_col(aes(y = pct), fill = "skyblue") +  # Bar plot for variance explained
  geom_line(aes(y = pct_cum, group = 1), color = "red") +  # Line plot for cumulative variance
  geom_point(aes(y = pct_cum), color = "red") +  # Points for cumulative variance
  labs(x = "Principal Component", y = "Fraction of Variance Explained (%)", 
       title = "Variance Explained by Principal Components") +
  theme_minimal()


#####################################################

# Load necessary libraries
library(dplyr)

# Extract loadings for PC1 and PC2
loadings_PC1 <- pca_result$rotation[, 1]  # PC1 loadings
loadings_PC2 <- pca_result$rotation[, 2]  # PC2 loadings

# Create data frames of loadings for PC1 and PC2
loadings_df_PC1 <- data.frame(Gene_ID = rownames(pca_result$rotation), PC1_Loading = loadings_PC1)
loadings_df_PC2 <- data.frame(Gene_ID = rownames(pca_result$rotation), PC2_Loading = loadings_PC2)

# Merge PC1 and PC2 loadings
loadings_df <- merge(loadings_df_PC1, loadings_df_PC2, by = "Gene_ID")

# Merge with gene_info to get gene names
loadings_df <- merge(loadings_df, gene_info, by = "Gene_ID", all.x = TRUE)

# Add a column to indicate whether the gene's contribution is positive or negative for PC1 and PC2
loadings_df <- loadings_df %>%
  mutate(PC1_Direction = ifelse(PC1_Loading > 0, "Positive", "Negative"),
         PC2_Direction = ifelse(PC2_Loading > 0, "Positive", "Negative"))

# Sort by absolute values of PC1 and PC2 loadings
loadings_df <- loadings_df %>%
  arrange(desc(abs(PC1_Loading)), desc(abs(PC2_Loading)))

# Save the list of genes to a CSV file
write.csv(loadings_df, "PC1_PC2_genes.csv", row.names = FALSE)

# Display the top genes for PC1 and PC2
print(loadings_df)


# Set a threshold to keep the top X% of genes based on absolute loadings for PC1 and PC2
threshold_percentage <- 0.10  # For example, top 10%

# Calculate the number of genes to keep
num_genes <- round(nrow(loadings_df) * threshold_percentage)

# Select the top genes based on absolute PC1 and PC2 loadings
top_genes <- loadings_df %>%
  top_n(num_genes, wt = abs(PC1_Loading)) %>%
  top_n(num_genes, wt = abs(PC2_Loading))

# Save the top contributing genes to a CSV file
write.csv(top_genes, "PC1_PC2_top10percent.csv", row.names = FALSE)

# Print the top genes
print(top_genes)

# Load necessary libraries
library(dplyr)
library(pheatmap)


# Extract expression values for these top 10% genes
expression_matrix <- dge_filtered$counts[top_genes$Gene_ID, ]

# Ensure rownames of the expression matrix are gene names for better labeling
rownames(expression_matrix) <- top_genes$Gene_Name

# Create distinct colors for each condition in the metadata using RColorBrewer
condition_colors <- brewer.pal(n = length(unique(metadata$Condition)), "Set1")
names(condition_colors) <- unique(metadata$Condition)

# Create an annotation dataframe for the conditions
annotation_col <- data.frame(Condition = metadata$Condition)
rownames(annotation_col) <- rownames(metadata)

# Create the heatmap with annotation and custom colors for the top 10% contributing genes
pheatmap(expression_matrix, 
         scale = "row",                         # Scale genes by row (Z-score)
         annotation_col = annotation_col,       # Annotate samples by condition
         annotation_colors = list(Condition = condition_colors), 
         show_rownames = FALSE,                  # Show gene names on heatmap
         main = "Heatmap of Top 10% Contributing Genes by Condition", 
         cluster_rows = FALSE,                   # Cluster rows (genes)
         cluster_cols = TRUE)                   # Cluster columns (samples)



######################################################

# Load necessary libraries
library(clusterProfiler)
library(org.Mm.eg.db)  # Mouse gene annotation package

# Assuming you have the Gene_Name column in your data frame
gene_list <- top_genes$Gene_Name

# Convert gene symbols to Entrez IDs if needed
entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Merge the Entrez IDs back into your original dataframe
top_contributing_genes_with_entrez <- merge(top_contributing_genes, entrez_ids, by.x = "Gene_Name", by.y = "SYMBOL", all.x = TRUE)

# Save the updated dataframe with Entrez IDs to a CSV file
write.csv(top_contributing_genes_with_entrez, "top_contributing_genes_with_EntrezIDs.csv", row.names = FALSE)

# Perform GO enrichment analysis
go_results <- enrichGO(gene         = entrez_ids$ENTREZID,
                       OrgDb        = org.Mm.eg.db,
                       keyType      = "ENTREZID",
                       ont          = "ALL",  # Can be "BP", "MF", "CC" or "ALL"
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

# Perform KEGG pathway enrichment analysis
kegg_results <- enrichKEGG(gene         = entrez_ids$ENTREZID,
                           organism     = 'mmu',  # 'mmu' is for mouse (Mus musculus)
                           pvalueCutoff = 0.05)

# Visualize top GO terms
dotplot(go_results, showCategory = 20) + 
  ggtitle("Top Enriched GO Terms for Genes Contributing to PC1 and PC2") +  # Only one ggtitle()
  theme(axis.text.y = element_text(size = 7))  # Adjust y-axis text size here

# Visualize top KEGG pathways
dotplot(kegg_results, showCategory = 20) + ggtitle("Top Enriched KEGG Pathways Associated with PC1 and PC2")+  # Only one ggtitle()
  theme(axis.text.y = element_text(size = 6)) 


###################################################

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

# Check the design matrix
head(design)

# Estimate the common dispersion
dge_filtered <- estimateGLMCommonDisp(dge_filtered, design)

# Estimate the trended dispersion
dge_filtered <- estimateGLMTrendedDisp(dge_filtered, design)

# Estimate the tagwise dispersion
dge_filtered <- estimateGLMTagwiseDisp(dge_filtered, design)

# Plot the biological coefficient of variation (BCV)
plotBCV(dge_filtered)


# Define contrasts for the comparisons you are interested in
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




# Loop through each contrast and apply the contrasts to the fitted model
for (contrast_name in names(contrasts_list)) {
  # Get the contrast matrix
  contrast <- contrasts_list[[contrast_name]]
  
  # Apply the contrasts to the fitted model
  fit_contrast <- contrasts.fit(fit, contrast)
  
  # Perform empirical Bayes moderation
  fit_contrast <- eBayes(fit_contrast)
  
  # Extract top differentially expressed genes for the contrast
  top_genes <- topTable(fit_contrast, number = Inf)
  
  # Save to a CSV file
  output_filename <- paste0(contrast_name, "_DEG_results.csv")
  write.csv(top_genes, file = output_filename, row.names = TRUE)
  
  # Print message to indicate saving progress
  cat("Saved results for", contrast_name, "to", output_filename, "\n")
  
  # Add a column for significance
  top_genes$Significant <- ifelse(top_genes$adj.P.Val < 0.05 & abs(top_genes$logFC) > 1, "Yes", "No")
  
  # Create the volcano plot
  p <- ggplot(top_genes, aes(x = logFC, y = -log10(P.Value), color = Significant)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("grey", "red")) +
    labs(title = paste0("Volcano Plot: ", contrast_name),
         x = "Log Fold Change",
         y = "-log10(P-Value)") +
    theme_minimal()
  
  # Print the plot
  print(p)
}

# Gene Ontology Enrichment Analysis
library(clusterProfiler)
library(org.Mm.eg.db)

# Loop through each contrast
for (contrast_name in names(contrasts_list)) {
  # Extract top differentially expressed genes for the contrast
  top_genes <- topTable(fit_contrast, number = Inf)
  
  # Filter significant genes based on adjusted p-value and log fold change
  significant_genes <- top_genes[top_genes$adj.P.Val < 0.05 & abs(top_genes$logFC) > 1, ]
  
  # Create a named vector for log fold changes
  gene_list <- significant_genes$logFC
  names(gene_list) <- rownames(significant_genes)
  
  # Perform GO enrichment analysis
  ego <- enrichGO(gene = names(gene_list), 
                  OrgDb = org.Mm.eg.db, 
                  keyType = "ENSEMBL",
                  ont = "BP",  # Biological Process
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.1)
  
  # Check if there are any significant GO terms
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    # Plot the results
    plot <- dotplot(ego, showCategory = 10) +
      ggtitle(paste0("GO Enrichment for ", contrast_name)) +
      theme(axis.text.y = element_text(size = 5))
    
    # Display the plot
    print(plot)
  } else {
    cat("No significant GO terms found for contrast:", contrast_name, "\n")
  }
}


# KEGG Pathway Analysis
library(ReactomePA)

# Convert the gene list to Entrez IDs (KEGG requires Entrez IDs)
entrez_ids <- bitr(condition_specific_genes, fromType = "ENSEMBL", 
                   toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Perform KEGG pathway enrichment analysis
kegg_enrich <- enrichKEGG(gene = entrez_ids$ENTREZID, 
                          organism = 'mmu',  # 'mmu' for mouse
                          pvalueCutoff = 0.05)

# Plot the KEGG results
dotplot(kegg_enrich, showCategory = 10) +
  theme(axis.text.y = element_text(size = 5))
