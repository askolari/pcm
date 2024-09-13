getwd() #should be pcm
setwd("C:/Users/agsko/dev/pcm/EdgeR")


library(here)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(limma)
library(edgeR)
library(openxlsx)
library(tibble)
library(pheatmap)
library(matrixStats)
library(SummarizedExperiment)
library(S4Vectors)

set.seed(1986)

# Data Loading and Preprocessing

## Loading Data

# Set the file paths to your local files
tpm_file <- "C:/Users/agsko/dev/pcm/TPM_values.csv"
raw_counts_file <- "C:/Users/agsko/dev/pcm/raw_counts.csv"
coldata_file <- "C:/Users/agsko/dev/pcm/coldata.xlsx"

# Load the CSV and Excel files
tpm_data <- read.csv(tpm_file)
raw_counts_data <- read.csv(raw_counts_file)
coldata <- read.xlsx(coldata_file)

# Display the first few rows of each dataset
head(tpm_data)
head(raw_counts_data)
head(coldata)

# Rename the column 'Sample' to 'Sample ID'
colnames(coldata)[colnames(coldata) == "Sample"] <- "Sample ID"

# Save the modified coldata back to the file
write.xlsx(coldata, "C:/Users/agsko/dev/pcm/coldata_modified.xlsx", rowNames = FALSE)

#Load the coldata_modified 
coldata_mod <- read.xlsx("C:/Users/agsko/dev/pcm/coldata_modified.xlsx")

# Load the additional metadata
sequencing_stats <- read.csv("C:/Users/agsko/dev/pcm/sequencing_stats.csv")
alignment_stats <- read.csv("C:/Users/agsko/dev/pcm/alignment_stats.csv")

# Combine the additional metadata with the coldata
# Assuming that the sequencing_stats and alignment_stats have matching sample IDs
metadata_combined <- merge(coldata_mod, sequencing_stats, by = "Sample.ID")
metadata_combined <- merge(metadata_combined, alignment_stats, by = "Sample.ID")

# Save the combined metadata back to the file
write.xlsx(metadata_combined, "C:/Users/agsko/dev/pcm/metadata_combined.xlsx", rowNames = FALSE)

# Display the structure of your datasets
str(raw_counts_data)
str(tpm_data)
str(metadata_combined)

# Extract gene names from the raw counts data
gene_names <- raw_counts_data$Gene.name


# Create a DataFrame to hold the gene names
gene_metadata <- DataFrame(Gene.name = gene_names)

# Reorder the columns: Move 'Gene.name' to the second position
raw_counts_data <- raw_counts_data %>% select(ID, Gene.name, everything())
tpm_data <- tpm_data %>% select(ID, Gene.name, everything())

# Create a rowData DataFrame from the 'ID' and 'Gene.name' columns
row_data <- DataFrame(ID = raw_counts_data$ID, Gene.name = raw_counts_data$Gene.name)

# Extract the counts and TPM matrices
counts_matrix <- as.matrix(raw_counts_data[, -c(1, 2)])# Exclude 'ID' and 'Gene.name'
rownames(counts_matrix) <- raw_counts_data$ID  # Set row names to gene IDs

# Check the first few row names of counts_matrix
head(rownames(counts_matrix))

tpm_matrix <- as.matrix(tpm_data[, -c(1, 2)])     # Exclude 'ID' and 'Gene.name'
rownames(tpm_matrix) <- tpm_data$ID  # Set row names to gene IDs

# Ensure that the sample names in metadata match the column names in counts_matrix and tpm_matrix
metadata_combined <- metadata_combined[match(colnames(counts_matrix), metadata_combined$Sample.ID), ]

# Check for NA values after matching
sum(is.na(metadata_combined$Sample.ID))

# Check the dimensions of row_data and counts_matrix
dim(row_data)
dim(counts_matrix)

# Check if the row names of counts_matrix match the ID column of row_data
identical(rownames(counts_matrix), row_data$ID)

# Check for duplicates in row_data$ID
sum(duplicated(row_data$ID))

# Check for duplicates in rownames(counts_matrix)
sum(duplicated(rownames(counts_matrix)))

# Check the result of the match() function
match_indices <- match(row_data$ID, rownames(counts_matrix))

# Check for any NA values in the match indices
sum(is.na(match_indices))

# Display the first few IDs from row_data and rownames(counts_matrix) for comparison
head(row_data$ID)
head(rownames(counts_matrix))

# Compare a few specific entries to see if there's a visible difference
row_data$ID[1:10]
rownames(counts_matrix)[1:10]

# Create the SummarizedExperiment object
sexp <- SummarizedExperiment(
  assays = list(counts = counts_matrix, TPM = tpm_matrix),
  rowData = row_data,
  colData = metadata_combined
)

# Save the SummarizedExperiment object if needed
saveRDS(sexp, "SummarizedExperiment_rnaseq_with_metadata.RDS")

# Load the SummarizedExperiment object from the .RDS file
sexp <- readRDS("SummarizedExperiment_rnaseq_with_metadata.RDS")

# Extract gene metadata from the SummarizedExperiment object
all.gene.meta <- as.data.frame(rowData(sexp))  # Extract rowData (gene metadata)

# Create a DGEList object from the raw counts
dge_test <- DGEList(counts = assay(sexp, "counts"))


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

#########################
#Variance partitioning

# Match samples based on 'Sample.ID' from metadata_combined and column names in dge_filtered
metadata_combined <- metadata_combined[match(colnames(dge_filtered), metadata_combined$Sample.ID), ]

# Check if sample IDs match between dge_filtered and metadata_combined
stopifnot(all(colnames(dge_filtered) == metadata_combined$Sample.ID))

# Voom transformation for linear modeling
vobj <- voom(dge_filtered, plot = FALSE)

# Define your formula for variance partitioning including sequencing quality and condition
# Modify this formula based on the column names in your metadata_combined
# e.g., Condition is your experimental condition, and Unique.Mapped.Reads.Level is your sequencing quality metric

form <- ~ Condition + Unique.Mapped.Reads + X..Bases....30  # Add/remove variables as needed

BiocManager::install("variancePartition", force = TRUE)
library(variancePartition)


# Check column names of expression matrix
colnames(vobj)

colnames(metadata_combined)


# Check row names of metadata
rownames(metadata_combined)

# Set row names of metadata_combined to the values in the Sample.ID column
rownames(metadata_combined) <- metadata_combined$Sample.ID

# Now check if the row names match the column names of vobj
identical(colnames(vobj), rownames(metadata_combined))  # This should return TRUE


# Perform variance partitioning
varPart <- fitExtractVarPartModel(vobj, form, metadata_combined)

# Sort variance partitioning results for visualization
vp <- sortCols(varPart)


# Generate the variance partitioning plot with smaller text
plotVarPart(vp, label.angle = 45) +
  ggtitle("Variance Partitioning by Condition and Sequencing Quality") +  # Add a proper title
  theme(
    axis.text.x = element_text(size = 8),  # Adjust x-axis text size
    axis.text.y = element_text(size = 8),  # Adjust y-axis text size
    axis.title.x = element_text(size = 10), # Adjust x-axis title size
    axis.title.y = element_text(size = 10), # Adjust y-axis title size
    plot.title = element_text(size = 14, hjust = 0.5)   # Adjust plot title size and center it
  )

# Save the variance partitioning plot
ggsave("var_partition_condition_quality.png")


#######################

#Data exploration based on multidimensional scaling
# Assuming 'metadata' contains the sample conditions
# Create a color vector based on the condition in your metadata

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
library(ggrepel)

# MDS plot for DGE after filtering and normalization
# Get the MDS coordinates from the dge_filtered object
mds_data <- plotMDS(dge_filtered, plot = FALSE)


# Convert to data frame
mds_coords <- data.frame(Sample = rownames(dge_filtered$samples),
                         x = mds_data$x, 
                         y = mds_data$y,
                         Condition = sample_conditions)  # Add the condition column

library(viridis)
library(RColorBrewer)
colors_set1 <- brewer.pal(n = length(unique(mds_coords$Condition)), name = "Set1")

# Plot using ggplot and ggrepel for label adjustment
ggplot(mds_coords, aes(x = x, y = y, label = Sample, color = Condition)) +
  geom_point(size = 3) +
  geom_text_repel(size = 2, nudge_y = 0.05, nudge_x = 0.05) +  # Adding nudge to avoid overlap
  theme_minimal() +
  labs(title = "MDS Plot After Normalisation", 
       x = paste0("Leading logFC dim 1 (", round(mds_data$var.explained[1] * 100, 1), "%)"), 
       y = paste0("Leading logFC dim 2 (", round(mds_data$var.explained[2] * 100, 1), "%)"),
       color = "Condition") +  # Adds label to the legend
  scale_color_manual(values = colors_set1) +  
  theme(legend.position = "right")

# Get the MDS coordinates from the dge_filtered object
mds_data <- plotMDS(dge_filtered, plot = FALSE)

# Convert to data frame
mds_coords <- data.frame(Sample = rownames(dge_filtered$samples),
                         x = mds_data$x, 
                         y = mds_data$y,
                         Condition = sample_conditions)  # Add the condition column


# Get a palette of pastel colors (you can use "Pastel1", "Pastel2", etc.)
pastel_colors <- brewer.pal(n = length(unique(mds_coords$Condition)), name = "Set2")

# Plot using ggplot and ggrepel for label adjustment
ggplot(mds_coords, aes(x = x, y = y, label = Sample, color = Condition)) +
  geom_point(size = 3) +
  geom_text_repel(size = 2, nudge_y = 0.05, nudge_x = 0.05) +  # Adding nudge to avoid overlap
  theme_minimal() +
  labs(title = "MDS Plot After Filtering", 
       x = paste0("Leading logFC dim 1 (", round(mds_data$var.explained[1] * 100, 1), "%)"), 
       y = paste0("Leading logFC dim 2 (", round(mds_data$var.explained[2] * 100, 1), "%)"),
       color = "Condition") +  # Adds label to the legend
  scale_color_manual(values = pastel_colors) +
  theme(legend.position = "right")  # Positions the legend on the side



# Select the top 500 most variable genes for clustering
gene_variances <- apply(dge_filtered$counts, 1, var)
top_500_genes <- names(sort(gene_variances, decreasing = TRUE))[1:500]
top_50_genes <- names(sort(gene_variances, decreasing = TRUE))[1:50]

# Subset the filtered counts for the top 500 genes
top500 <- dge_filtered$counts[top_500_genes, ]
top50 <- dge_filtered$counts[top_50_genes, ]

# Create the heatmap with hierarchical clustering for the top 500 genes
pheatmap(top500,
         scale = "row",                           # Scale genes by row (Z-score)
         annotation_col = metadata_combined[c("Condition")],  # Annotate samples by condition
         annotation_colors = list(Condition = condition_colors), 
         show_rownames = FALSE,                   # Hide gene names on heatmap
         main = "Heatmap of Top 500 Most Variable Genes by Condition", 
         cluster_rows = TRUE,                     # Cluster rows (genes)
         cluster_cols = TRUE)                     # Cluster columns (samples)

# Create the heatmap with hierarchical clustering for the top 500 genes
pheatmap(top50,
         scale = "row",                           # Scale genes by row (Z-score)
         annotation_col = metadata_combined[c("Condition")],  # Annotate samples by condition
         annotation_colors = list(Condition = condition_colors), 
         show_rownames = TRUE,                   # Hide gene names on heatmap
         main = "Heatmap of Top 50 Most Variable Genes by Condition", 
         cluster_rows = TRUE,                     # Cluster rows (genes)
         cluster_cols = TRUE)  

#####


# Check variance in expression for all samples
expression_variance <- apply(dge_filtered$counts, 2, var)


# Print variance of X94 compared to X92 and X93
print(expression_variance[c("X92", "X93", "X94", "X97", "X99",  "X101")])

#X92      X93      X94      X97      X99     X101 
#26394813 24203102 18167424 69400596 15266400 31194829

# Perform PCA on normalized counts (after filtering)
pca_result <- prcomp(t(dge_filtered$counts), scale. = TRUE)

# Create a data frame for PCA results
pca_data <- as.data.frame(pca_result$x)
pca_data$Sample <- rownames(pca_data)

# Add the sample conditions to the PCA data
pca_data$Condition <- metadata$Condition

percent_var <- round(100 * summary(pca_result)$importance[2, 1:2], 2) # for PC1 and PC2

# Modify the plot
ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(aes(color=Condition)) +
  geom_text_repel(aes(label = Sample), show.legend = FALSE) +  # Use geom_text_repel
  theme_minimal() +
  ggtitle("PCA Plot of Samples by Condition") +
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

ggplot(pc_eigenvalues, aes(x = PC)) +
  geom_col(aes(y = pct), fill = "skyblue") +  # Bar plot for variance explained
  geom_line(aes(y = pct_cum, group = 1), color = "red") +  # Line plot for cumulative variance
  geom_point(aes(y = pct_cum), color = "red") +  # Points for cumulative variance
  labs(x = "Principal Component", y = "Fraction of Variance Explained (%)", 
       title = "Variance Explained by Principal Components") +
  theme_minimal()

# Install and load necessary packages
if (!requireNamespace("plotly", quietly = TRUE)) {
  install.packages("plotly")
}

library(plotly)

# Perform PCA on the filtered counts data
pca_result <- prcomp(t(dge_filtered$counts), scale. = TRUE)

# Extract the first three principal components
pca_data <- as.data.frame(pca_result$x[, 1:3])
colnames(pca_data) <- c("PC1", "PC2", "PC3")

# Add sample conditions to the PCA data for coloring
pca_data$Condition <- metadata_combined$Condition  # Assuming metadata_combined contains sample condition info

# Create a 3D PCA plot using plotly
p <- plot_ly(pca_data, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Condition, colors = c("blue", "red", "green", "purple", "orange"),
             type = "scatter3d", mode = "markers", marker = list(size = 6))

# Customize the layout of the 3D plot
p <- p %>%
  layout(title = "3D PCA Plot",
         scene = list(
           xaxis = list(title = "PC1"),
           yaxis = list(title = "PC2"),
           zaxis = list(title = "PC3")
         ))

# Show the plot
p



#####################################################

# Load necessary libraries

# Extract loadings for PC1 and PC2
loadings_PC1 <- pca_result$rotation[, 1]  # PC1 loadings
loadings_PC2 <- pca_result$rotation[, 2]  # PC2 loadings

# Create data frames of loadings for PC1 and PC2
loadings_df_PC1 <- data.frame(Gene_ID = rownames(pca_result$rotation), PC1_Loading = loadings_PC1)
loadings_df_PC2 <- data.frame(Gene_ID = rownames(pca_result$rotation), PC2_Loading = loadings_PC2)

# Merge PC1 and PC2 loadings
loadings_df <- merge(loadings_df_PC1, loadings_df_PC2, by = "Gene_ID")

# Merge with gene_info to get gene names
loadings_df <- merge(loadings_df, raw_counts, by = "Gene_ID", all.x = TRUE)

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
threshold_percentage <- 0.005 

# Calculate the number of genes to keep
num_genes <- round(nrow(loadings_df) * threshold_percentage)

# Select the top genes based on absolute PC1 and PC2 loadings
top_genes <- loadings_df %>%
  top_n(num_genes, wt = abs(PC1_Loading)) %>%
  top_n(num_genes, wt = abs(PC2_Loading))

# Save the top contributing genes to a CSV file
write.csv(top_genes, "PC1_PC2_p0.005.csv", row.names = FALSE)

# Print the top genes
print(top_genes)

# Check if rownames of metadata_combined are sample names
rownames(metadata_combined) <- metadata_combined$Sample.ID

# Define condition colors
condition_colors <- c(
  "Ure" = "blue",
  "Ure_DMXAA" = "red",
  "Ure_RT" = "green",
  "Ure_RT_DMXAA" = "purple",
  "PBS_RT_DMXAA" = "orange"
)


#Select the top 20 genes contributing to PC1 and PC2
top_20_genes <- loadings_df[order(abs(loadings_df$PC1_Loading) + abs(loadings_df$PC2_Loading), decreasing = TRUE)[1:20], "Gene_ID"]

#Map Ensembl IDs to gene names using 'all.gene.meta'
# 'all.gene.meta' has columns 'ID' (Ensembl) and 'Gene.name'
top_20_genes_with_names <- merge(data.frame(Gene_ID = top_20_genes), all.gene.meta, by.x = "Gene_ID", by.y = "ID")

#Subset the filtered counts for the top 20 genes based on Ensembl IDs
top20 <- dge_filtered$counts[top_20_genes_with_names$Gene_ID, ]
rownames(top20) <- top_20_genes_with_names$Gene.name  # Replace row names with gene names

# Step 4: Create the heatmap for the top 20 genes
pheatmap(top20,
         scale = "row",                           # Scale genes by row (Z-score)
         annotation_col = metadata_combined[c("Condition")],  # Annotate samples by condition
         annotation_colors = list(Condition = condition_colors), 
         show_rownames = TRUE,                    # Show gene names on heatmap
         main = "Heatmap of Top 20 Genes from PC1 and PC2 by Condition", 
         cluster_rows = TRUE,                     # Cluster rows (genes)
         cluster_cols = TRUE) 


##############

# Sort the loadings by absolute values of PC1 and PC2
loadings_df$PC1_Loading <- abs(loadings_df$PC1_Loading)
loadings_df$PC2_Loading <- abs(loadings_df$PC2_Loading)

top20_pc1 <- loadings_df %>%
  arrange(desc(PC1_Loading)) %>%
  slice(1:20) %>%
  dplyr::select(Gene_ID, PC1_Loading)

# Get the top 20 genes contributing to PC2
top20_pc2 <- loadings_df %>%
  arrange(desc(PC2_Loading)) %>%
  slice(1:20) %>%
  dplyr::select(Gene_ID, PC2_Loading)

# Merge the top 20 PC1 genes with gene names
top20_pc1_with_names <- merge(top20_pc1, raw_counts_data[, c("ID", "Gene.name")], 
                              by.x = "Gene_ID", by.y = "ID", all.x = TRUE)

# Merge the top 20 PC2 genes with gene names
top20_pc2_with_names <- merge(top20_pc2, raw_counts_data[, c("ID", "Gene.name")], 
                              by.x = "Gene_ID", by.y = "ID", all.x = TRUE)

# Print the top 20 PC1 and PC2 genes with names
cat("Top 20 genes contributing to PC1 with Gene Names:\n")
print(top20_pc1_with_names)

cat("\nTop 20 genes contributing to PC2 with Gene Names:\n")
print(top20_pc2_with_names)

# Save the top 20 genes from PC1
write.csv(top20_pc1_with_names, "top20_PC1_genes_with_names.csv", row.names = FALSE)

# Save the top 20 genes from PC2
write.csv(top20_pc2_with_names, "top20_PC2_genes_with_names.csv", row.names = FALSE)

# Combine the top 20 genes from PC1 and PC2
top_genes_pc1_pc2 <- unique(c(top20_pc1_with_names$Gene_ID, top20_pc2_with_names$Gene_ID))


######################################################

gene_list <- all.gene.meta$Gene.name

# Merge top_genes with all_gene_meta based on the gene IDs
top_genes <- merge(top_genes, all.gene.meta, by.x = "Gene_ID", by.y = "ID", all.x = TRUE)

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Merge the Entrez IDs back into original df
top_contributing_genes_with_entrez <- merge(top_genes, entrez_ids, by.x = "Gene.name.x", by.y = "SYMBOL", all.x = TRUE)

# Save the updated dataframe with Entrez IDs to a CSV file
write.csv(top_contributing_genes_with_entrez, "top_contributing_genes_with_EntrezIDs.csv", row.names = FALSE)

# Perform GO enrichment analysis
go_results <- enrichGO(gene         = entrez_ids$ENTREZID,
                       OrgDb        = org.Mm.eg.db,
                       keyType      = "ENTREZID",
                       ont          = "ALL",  # Can be "BP", "MF", "CC" or "ALL"
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.1)

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

#######################
#K-means

library(factoextra)  # For clustering visualization

# Assuming `voom_logcpm` contains the logCPM values from voom-transformation
# Or, if you're using normalized counts, use: dge_filtered$counts (normalized)

# Select your data for clustering (e.g., logCPM values)
clustering_data <- t(dge_filtered$counts)  # Transpose and Use your voom logCPM or normalized counts matrix

# Perform k-means clustering with k = 3 (start with 3, but you can iterate for best k)
set.seed(123)
kmeans_result <- kmeans(clustering_data, centers = 4, nstart = 25)

# Add cluster assignments to the metadata or a new dataframe
metadata_combined$Cluster <- as.factor(kmeans_result$cluster)

# Visualize K-means clustering using PCA or MDS for dimensionality reduction

# Run PCA on logCPM data for visualization
pca_result <- prcomp(clustering_data, scale. = TRUE)

# Create a dataframe with PCA results
pca_data <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2], Cluster = metadata_combined$Cluster)

# Plot the clusters using ggplot
ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 3) +
  labs(title = "K-means Clustering (k = 4) on PCA-transformed Data") +
  theme_minimal()

# Optionally, use a scree plot to select the optimal number of clusters (elbow method)

fviz_nbclust(clustering_data, kmeans, method = "wss") +
  labs(subtitle = "Elbow method for optimal k")

# Create a data frame that combines sample IDs, conditions, and cluster assignments
cluster_assignments <- data.frame(
  Sample_ID = rownames(metadata_combined),  # Use the sample IDs from metadata_combined
  Condition = metadata_combined$Condition,  # Use the condition information
  Cluster = kmeans_result$cluster  # The cluster assignments from k-means
)

# View the cluster assignments table
print(cluster_assignments)
# Save the cluster assignments to a CSV file
write.csv(cluster_assignments, "cluster_assignments.csv", row.names = FALSE)

###########
#DGE within clusters
# Subset the DGEList object for each cluster
# Subset the Condition for Cluster 1 and ensure it's a factor
metadata_combined$Condition <- as.factor(metadata_combined$Condition)
# Subset the original DGEList based on samples in Cluster 1
dge_cluster1 <- dge_filtered[, metadata_combined$Cluster == 1]

# Check if the samples are correctly subset
dim(dge_cluster1)


condition_cluster1 <- factor(metadata_combined$Condition[metadata_combined$Cluster == 1])


# Check the levels (the first level is the baseline)
levels(condition_cluster1)

# Relevel to set "Ure" as the baseline
condition_cluster1 <- relevel(condition_cluster1, ref = "Ure_RT")

# Create the design matrix with "Ure" as the baseline
design_cluster1 <- model.matrix(~ condition_cluster1)

# Estimate dispersion and fit the model
dge_cluster1 <- estimateDisp(dge_cluster1, design_cluster1)
fit_cluster1 <- glmFit(dge_cluster1, design_cluster1)
lrt_cluster1 <- glmLRT(fit_cluster1)

# Extract top differentially expressed genes
top_genes_cluster1 <- topTags(lrt_cluster1, n = 100)
print(top_genes_cluster1)

# Save the results
write.csv(lrt_cluster1, "dge_cluster1.csv")



# Subset metadata for Cluster 2
cluster_2 <- metadata_combined[metadata_combined$Cluster == 2, ]

# Subset dge_filtered for sample X97
dge_cluster_2 <- dge_filtered[, colnames(dge_filtered) == "X97"]

# Check if the subset has the correct count data
colnames(dge_cluster_2$counts)

# Extract gene IDs and gene names from dge_filtered
gene_ids <- rownames(dge_filtered$counts)
gene_names <- dge_filtered$genes$Gene.name

# Check if the gene names match the gene IDs
head(gene_ids)
head(gene_names)

# Create a data frame with gene IDs, gene names, and counts for X97
cluster_2_gene_data <- data.frame(
  Gene_ID = gene_ids,
  Gene_Name = gene_names,
  Count = dge_cluster_2$counts[, 1]  # Get the count for the single sample in Cluster 2 (X97)
)

# Check the first few rows of the data frame
head(cluster_2_gene_data)

# Save the gene data for cluster 2 to a CSV file
write.csv(cluster_2_gene_data, "cluster_2_genes.csv", row.names = FALSE)

##
# Subset the DGEList for Cluster 3 samples
cluster_3_samples <- c("X88", "X89", "X90", "X92", "X93", "X94")
dge_cluster_3 <- dge_filtered[, colnames(dge_filtered) %in% cluster_3_samples]

# Ensure the normalization factors are inherited from the full dataset
dge_cluster_3$samples$norm.factors <- dge_filtered$samples$norm.factors[colnames(dge_filtered) %in% cluster_3_samples]

# Check if normalization factors are properly transferred
print(dge_cluster_3$samples$norm.factors)

# Check normalization factors to ensure they are not NA
if (any(is.na(dge_cluster_3$samples$norm.factors))) {
  stop("Normalization factors are still missing.")
}

# Create a group variable for Cluster 3 samples
group_cluster_3 <- factor(c(rep("Ure", 3), rep("Ure_DMXAA", 3)))

# Estimate dispersion
dge_cluster_3 <- estimateDisp(dge_cluster_3, design = model.matrix(~ group_cluster_3))

# Fit the GLM model
fit_cluster_3 <- glmFit(dge_cluster_3, design = model.matrix(~ group_cluster_3))

# Perform the likelihood ratio test
lrt_cluster_3 <- glmLRT(fit_cluster_3, coef = 2)

condition_cluster3 <- factor(metadata_combined$Condition[metadata_combined$Cluster == 3])

# Check the levels (the first level is the baseline)
levels(condition_cluster3)

# Get the top differentially expressed genes
top_genes_cluster_3 <- topTags(lrt_cluster_3, n = 20)
print(top_genes_cluster_3)

# Extract all differential expression results within Cluster 3
all_genes_cluster_3 <- lrt_cluster_3$table

# Save the full results with gene names for Cluster 3
write.csv(all_genes_cluster_3, "dge_within_cluster3.csv", row.names = TRUE)


condition_cluster3 <- factor(metadata_combined$Condition[metadata_combined$Cluster == 3])

# Create the design matrix with "Ure" as the baseline
design_cluster3 <- model.matrix(~ condition_cluster3)

# Subset the DGEList object for Cluster 3 samples
cluster_3_samples <- colnames(dge_filtered)[metadata_combined$Cluster == 3]
dge_cluster3 <- dge_filtered[, cluster_3_samples]

# Ensure the normalization factors are inherited from the full dataset
dge_cluster3$samples$norm.factors <- dge_filtered$samples$norm.factors[colnames(dge_filtered) %in% cluster_3_samples]

# Check if normalization factors are properly transferred
print(dge_cluster3$samples$norm.factors)

# Estimate dispersion for Cluster 4
dge_cluster3 <- estimateDisp(dge_cluster3, design_cluster3)

# Check the dispersion estimates
print(dge_cluster3$common.dispersion)

# Fit the GLM model for Cluster 4
fit_cluster3 <- glmFit(dge_cluster3, design_cluster3)

# Perform the likelihood ratio test (LRT) for the comparison with the baseline "Ure_RT_DMXAA"
lrt_cluster3_ure_rt_dmxaa <- glmLRT(fit_cluster3)

# Extract top differentially expressed genes for "PBS_RT_DMXAA vs Ure_RT_DMXAA"
top_genes_cluster3_ure_vs_ure_dmxaa <- topTags(lrt_cluster3_ure_rt_dmxaa, n = Inf)$table

# Save the results
write.csv(top_genes_cluster3_ure_vs_ure_dmxaa, "dge_cluster3_ure_vs_ure_dmxaa.csv", row.names = TRUE)


###########################
#Upregulated and downregulated genes
# Load required libraries
library(dplyr)

# Define the file paths
file_paths <- c(
  "C:/Users/agsko/dev/pcm/EdgeR/dge_cluster3_ure_vs_ure_dmxaa.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/dge_cluster4_pbs_vs_ure_rt.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/dge_cluster4_ure_rt_vs_ure_rt_dmxaa.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/dge_cluster4_pbs_vs_ure_rt_dmxaa.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/dge_cluster1_ure_rt_vs_pbs_rt_dmxaa.csv"
)

# Function to filter upregulated and downregulated genes, using PValue when FDR is missing
filter_genes_with_pvalue <- function(file_path) {
  # Read the CSV file
  data <- read.csv(file_path)
  
  # Check if 'FDR' column exists
  if ("FDR" %in% colnames(data)) {
    # Filter upregulated genes (logFC > 1 & FDR < 0.001)
    upregulated_genes <- data[data$logFC > 1 & data$FDR < 0.001, ]
    
    # Filter downregulated genes (logFC < -1 & FDR < 0.001)
    downregulated_genes <- data[data$logFC < -1 & data$FDR < 0.001, ]
  } else {
    # Filter using PValue if FDR is missing (logFC > 1 & PValue < 0.001)
    upregulated_genes <- data[data$logFC > 1 & data$PValue < 0.001, ]
    
    # Filter downregulated genes (logFC < -1 & PValue < 0.001)
    downregulated_genes <- data[data$logFC < -1 & data$PValue < 0.001, ]
  }
  
  # Save the filtered genes to new CSV files
  file_base <- gsub(".csv", "", basename(file_path)) # Get the file base name
  up_file <- paste0(file_base, "_upregulated.csv")
  down_file <- paste0(file_base, "_downregulated.csv")
  
  write.csv(upregulated_genes, paste0("C:/Users/agsko/dev/pcm/EdgeR/", up_file), row.names = FALSE)
  write.csv(downregulated_genes, paste0("C:/Users/agsko/dev/pcm/EdgeR/", down_file), row.names = FALSE)
  
  return(list(upregulated_genes = upregulated_genes, downregulated_genes = downregulated_genes))
}

# Apply the function to all files
results_with_pvalue <- lapply(file_paths, filter_genes_with_pvalue)




metadata_combined$Condition <- as.factor(metadata_combined$Condition)

condition_cluster4 <- factor(metadata_combined$Condition[metadata_combined$Cluster == 4])


# Check the levels (the first level is the baseline)
levels(condition_cluster4)

# Relevel to set "Ure" as the baseline
condition_cluster4 <- relevel(condition_cluster4, ref = "Ure_RT_DMXAA")

# Create the design matrix with "Ure" as the baseline
design_cluster4 <- model.matrix(~ condition_cluster4)

# Subset the DGEList object for Cluster 4 samples
cluster_4_samples <- colnames(dge_filtered)[metadata_combined$Cluster == 4]
dge_cluster4 <- dge_filtered[, cluster_4_samples]

# Ensure the normalization factors are inherited from the full dataset
dge_cluster4$samples$norm.factors <- dge_filtered$samples$norm.factors[colnames(dge_filtered) %in% cluster_4_samples]

# Check if normalization factors are properly transferred
print(dge_cluster4$samples$norm.factors)

# Estimate dispersion for Cluster 4
dge_cluster4 <- estimateDisp(dge_cluster4, design_cluster4)

# Check the dispersion estimates
print(dge_cluster4$common.dispersion)

# Fit the GLM model for Cluster 4
fit_cluster4 <- glmFit(dge_cluster4, design_cluster4)

# Perform the likelihood ratio test (LRT) for the comparison with the baseline "Ure_RT_DMXAA"
lrt_cluster4_ure_rt_dmxaa <- glmLRT(fit_cluster4, coef = 2)

# Extract top differentially expressed genes for "PBS_RT_DMXAA vs Ure_RT_DMXAA"
top_genes_cluster4_ure_rt_dmxaa <- topTags(lrt_cluster4_ure_rt_dmxaa, n = Inf)$table

# Save the results
write.csv(top_genes_cluster4_ure_rt_dmxaa, "dge_cluster4_ure_rt_dmxaa.csv", row.names = TRUE)

# Coefficients for the design matrix
# coef = 2 compares "PBS_RT_DMXAA" to "Ure_RT_DMXAA"
# coef = 3 compares "Ure_RT" to "Ure_RT_DMXAA"

# Perform the LRT for "PBS_RT_DMXAA" vs "Ure_RT_DMXAA"
lrt_cluster4_pbs_vs_ure_rt_dmxaa <- glmLRT(fit_cluster4, coef = 2)

# Perform the LRT for "Ure_RT" vs "Ure_RT_DMXAA"
lrt_cluster4_ure_rt_vs_ure_rt_dmxaa <- glmLRT(fit_cluster4, coef = 3)

# Extract top differentially expressed genes for both comparisons
top_genes_cluster4_pbs_vs_ure_rt_dmxaa <- topTags(lrt_cluster4_pbs_vs_ure_rt_dmxaa, n = Inf)$table
top_genes_cluster4_ure_rt_vs_ure_rt_dmxaa <- topTags(lrt_cluster4_ure_rt_vs_ure_rt_dmxaa, n = Inf)$table


# Save the results for "PBS_RT_DMXAA vs Ure_RT_DMXAA"
write.csv(top_genes_cluster4_pbs_vs_ure_rt_dmxaa, "dge_cluster4_pbs_vs_ure_rt_dmxaa.csv", row.names = TRUE)

# Save the results for "Ure_RT vs Ure_RT_DMXAA"
write.csv(top_genes_cluster4_ure_rt_vs_ure_rt_dmxaa, "dge_cluster4_ure_rt_vs_ure_rt_dmxaa.csv", row.names = TRUE)

# Relevel to set "Ure_RT" as the baseline
condition_cluster4 <- relevel(condition_cluster4, ref = "Ure_RT")

# Create the design matrix with "Ure_RT" as the baseline
design_cluster4 <- model.matrix(~ condition_cluster4)

# Fit the GLM model for Cluster 4 with the new baseline
dge_cluster4 <- estimateDisp(dge_cluster4, design_cluster4)
fit_cluster4 <- glmFit(dge_cluster4, design_cluster4)

# Perform the likelihood ratio test (LRT) for the comparison between "PBS_RT_DMXAA" and "Ure_RT"
lrt_cluster4_pbs_vs_ure_rt <- glmLRT(fit_cluster4, coef = 2)

# Extract top differentially expressed genes for "PBS_RT_DMXAA vs Ure_RT"
top_genes_cluster4_pbs_vs_ure_rt <- topTags(lrt_cluster4_pbs_vs_ure_rt, n = Inf)$table

# Save the results to a CSV file
write.csv(top_genes_cluster4_pbs_vs_ure_rt, "dge_cluster4_pbs_vs_ure_rt.csv", row.names = TRUE)


###################################################
#enter entrez ID


library(biomaRt)
library(dplyr)

# Set up Biomart for mouse genes
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Function to add Entrez IDs to a file based on Gene.name
add_entrez_id <- function(file_path) {
  # Read the data
  data <- read.csv(file_path)
  
  # Check for missing values in Gene.name and remove them
  data <- data %>% filter(!is.na(Gene.name) & Gene.name != "")
  
  # Retrieve Entrez IDs using Gene.name
  gene_names <- unique(data$Gene.name)
  gene_mapping <- getBM(
    attributes = c("mgi_symbol", "entrezgene_id"),
    filters = "mgi_symbol",
    values = gene_names,
    mart = ensembl
  )
  
  # Check for many-to-many relationships
  if (nrow(gene_mapping) != length(gene_names)) {
    warning("Detected a potential many-to-many relationship in gene mapping")
  }
  
  # Merge with the original data, handling potential duplicates
  merged_data <- left_join(data, gene_mapping, by = c("Gene.name" = "mgi_symbol"), relationship = "many-to-many")
  
  # Save the updated file with Entrez IDs
  new_file_path <- gsub(".csv", "_with_entrez.csv", file_path)
  write.csv(merged_data, new_file_path, row.names = FALSE)
  
  return(new_file_path)
}

# List of file paths to process
file_paths <- c(
  "C:/Users/agsko/dev/pcm/EdgeR/dge_cluster1_ure_rt_vs_pbs_rt_dmxaa.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/cluster_2_genes.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/dge_cluster1_ure_rt_vs_pbs_rt_dmxaa_downregulated.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/dge_cluster1_ure_rt_vs_pbs_rt_dmxaa_upregulated.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/dge_cluster4_pbs_vs_ure_rt_dmxaa_downregulated.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/dge_cluster4_pbs_vs_ure_rt_dmxaa_upregulated.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/dge_cluster4_ure_rt_vs_ure_rt_dmxaa_downregulated.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/dge_cluster4_ure_rt_vs_ure_rt_dmxaa_upregulated.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/dge_cluster4_pbs_vs_ure_rt_downregulated.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/dge_cluster4_pbs_vs_ure_rt_upregulated.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/dge_cluster3_ure_vs_ure_dmxaa_downregulated.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/dge_cluster3_ure_vs_ure_dmxaa_upregulated.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/dge_cluster3_ure_vs_ure_dmxaa.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/dge_cluster4_pbs_vs_ure_rt.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/dge_cluster4_ure_rt_vs_ure_rt_dmxaa.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/dge_cluster4_pbs_vs_ure_rt_dmxaa.csv"
)

# Apply the function to all files
lapply(file_paths, add_entrez_id)

########################
#GO and FGSEA


# Load all pathway files
mh_pathways <- gmtPathways("C:/Users/agsko/dev/pcm/mh.all.v2024.1.Mm.entrez.gmt")
m8_pathways <- gmtPathways("C:/Users/agsko/dev/pcm/m8.all.v2024.1.Mm.entrez.gmt")
m5_pathways <- gmtPathways("C:/Users/agsko/dev/pcm/m5.all.v2024.1.Mm.entrez.gmt")
m3_pathways <- gmtPathways("C:/Users/agsko/dev/pcm/m3.all.v2024.1.Mm.entrez.gmt")
m2_pathways <- gmtPathways("C:/Users/agsko/dev/pcm/m2.all.v2024.1.Mm.entrez.gmt")


# Read in the DGE data
dge_file <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/dge_cluster3_ure_vs_ure_dmxaa.csv")

levels(dge_file)

# Remove rows with NA in the 'entrezgene_id' column
dge_file <- na.omit(dge_file, cols = "entrezgene_id")

# Make sure the gene list is sorted by logFC
dge_file <- dge_file[order(-dge_file$logFC), ]

# Create a named vector of logFC, where the names are the Entrez IDs
ranks <- setNames(dge_file$logFC, dge_file$entrezgene_id)

library(fgsea)

# Perform FGSEA
fgsea_res <- fgsea(pathways = mh_pathways, 
                   stats = ranks, 
                   minSize = 15, 
                   maxSize = 500 
                   )

# Filter results by adjusted p-value
fgsea_res <- fgsea_res[fgsea_res$padj < 0.05, ]

# Plot the results
fgsea_plot <- ggplot(fgsea_res, aes(reorder(pathway, NES), NES)) +
  geom_bar(stat = "identity", aes(fill = padj), show.legend = TRUE) +
  scale_fill_gradient(low = "blue", high = "red") +
  coord_flip() +
  labs(title = "Enriched Pathways (GSEA) for Cluster 4(MH)\n(Ure_RT vs Ure_RT_DMXAA)", 
       x = "Pathway", 
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal()

# Show the plot
print(fgsea_plot)

# Save the plot if needed
ggsave("fgsea_cluster1_m8_pathways.png", fgsea_plot)




############

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
  Ure_DMXAA_vs_Ure = makeContrasts(Ure_DMXAA - Ure, levels = design),
  Ure_RT_vs_Ure = makeContrasts(Ure_RT - Ure, levels = design),
  Ure_RT_DMXAA_vs_Ure_RT = makeContrasts(Ure_RT_DMXAA - Ure_RT, levels = design),
  Ure_RT_DMXAA_vs_Ure = makeContrasts(Ure_RT_DMXAA - Ure, levels = design),
  Ure_vs_PBS_RT_DMXAA = makeContrasts(Ure - PBS_RT_DMXAA, levels = design),
  Ure_DMXAA_vs_PBS_RT_DMXAA = makeContrasts(Ure_DMXAA - PBS_RT_DMXAA, levels = design),
  Ure_RT_vs_PBS_RT_DMXAA = makeContrasts(Ure_RT - PBS_RT_DMXAA, levels = design),
  Ure_RT_DMXAA_vs_PBS_RT_DMXAA = makeContrasts(Ure_RT_DMXAA - PBS_RT_DMXAA, levels = design),
  Ure_DMXAA_vs_Ure_RT = makeContrasts(Ure_DMXAA - Ure_RT, levels = design),
  Ure_RT_DMXAA_vs_Ure_DMXAA = makeContrasts(Ure_RT_DMXAA - Ure_DMXAA, levels = design)
)

####

# Fit the GLM model to the data
fit <- glmFit(dge_filtered, design)

# Loop through each contrast and save the results
for (contrast_name in names(contrasts_list)) {
  # Get the contrast matrix
  contrast <- contrasts_list[[contrast_name]]
  
  # Apply the contrast and perform the LRT
  lrt <- glmLRT(fit, contrast = contrast)
  
  # Step 1: Extract the top 15,000 differentially expressed genes based on p-value ranking
  top_genes_all <- topTags(lrt, n = 15000)$table  # Extract top 15,000 genes based on p-value ranking
  
  # Save the result for all genes (up to 15,000)
  output_filename_all <- paste0(contrast_name, "_DEG_15,000_genes.csv")
  write.csv(top_genes_all, file = output_filename_all, row.names = TRUE)
  
  # Step 2: Apply a p-value cutoff (e.g., p < 0.001) to filter for significant genes
  top_genes_significant <- top_genes_all[top_genes_all$PValue < 0.001, ]  # Filter by p-value
  
  # Save the result for significant genes
  output_filename_sig <- paste0(contrast_name, "_DEG_p0.0001_genes.csv")
  write.csv(top_genes_significant, file = output_filename_sig, row.names = TRUE)
  
  # Print message to indicate saving progress
  cat("Saved results for", contrast_name, "\n")
}

# Loop through each contrast and save the results
for (contrast_name in names(contrasts_list)) {
  # Get the contrast matrix
  contrast <- contrasts_list[[contrast_name]]
  
  # Apply the contrast and perform the LRT
  lrt <- glmLRT(fit, contrast = contrast)
  
  # Step 1: Extract the top 15,000 differentially expressed genes based on p-value ranking
  top_genes_all <- topTags(lrt, n = 20000)$table  # Extract top 15,000 genes based on p-value ranking
  
  # Save the result for all genes (up to 15,000)
  output_filename_all <- paste0(contrast_name, "_DEG_20,000_genes.csv")
  write.csv(top_genes_all, file = output_filename_all, row.names = TRUE)
  
  # Print message to indicate saving progress
  cat("Saved results for", contrast_name, "\n")
}
###########################

# Loop through each contrast and save the results
for (contrast_name in names(contrasts_list)) {
  # Get the contrast matrix
  contrast <- contrasts_list[[contrast_name]]
  
  # Apply the contrast and perform the LRT
  lrt <- glmLRT(fit, contrast = contrast)
  
  # Step 3: Use decideTestsDGE to identify DEGs
  deGenes <- decideTestsDGE(lrt, p.value = 0.001)
  
  # Extract the gene names of significant DEGs
  deGenes_list <- rownames(lrt)[as.logical(deGenes)]
  
  # Step 4: Create a smear plot to visualize DEGs
  plotSmear(lrt, de.tags = deGenes_list)
  abline(h = c(-1, 1), col = 2)  # Add a horizontal line at log-fold change ±1
  
  # Add title to the plot with contrast name and p-value threshold
  title(main = paste(contrast_name, "- Smear Plot", "\n", "p < 0.001"))
  
  # Optionally, save the plot as an image
  plot_filename <- paste0(contrast_name, "_smear_plot.png")
  png(plot_filename)
  plotSmear(lrt, de.tags = deGenes_list)
  abline(h = c(-1, 1), col = 2)
  title(main = paste(contrast_name, "- Smear Plot", "\n", "p < 0.001"))
  dev.off()
  
  # Print message indicating the plot was saved
  cat("Saved smear plot for", contrast_name, "\n")
}

library(ggplot2)

# Loop through each contrast and save the results
for (contrast_name in names(contrasts_list)) {
  # Get the contrast matrix
  contrast <- contrasts_list[[contrast_name]]
  
  # Apply the contrast and perform the LRT
  lrt <- glmLRT(fit, contrast = contrast)
  
  # Step 3: Use decideTestsDGE to identify DEGs
  deGenes <- decideTestsDGE(lrt, p.value = 0.001)
  
  # Extract the top genes from the LRT result
  top_genes <- topTags(lrt, n = Inf)$table
  
  # Add a column to indicate significance based on p-value and log-fold change
  top_genes$Significant <- ifelse(top_genes$PValue < 0.001 & abs(top_genes$logFC) > 1, "Yes", "No")
  
  # Step 4: Create a volcano plot
  p <- ggplot(top_genes, aes(x = logFC, y = -log10(PValue), color = Significant)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("grey", "red")) +
    labs(title = paste(contrast_name, "- Volcano Plot", "\n", "p < 0.001"),
         x = "Log Fold Change",
         y = "-log10(P-Value)") +
    theme_minimal()
  
  # Display the plot
  print(p)
  
  # Optionally, save the plot as an image
  plot_filename <- paste0(contrast_name, "_volcano_plot.png")
  ggsave(plot_filename, plot = p, width = 6, height = 6)
  
  # Print message indicating the plot was saved
  cat("Saved volcano plot for", contrast_name, "\n")
}

################################

# Loop through each contrast
for (contrast_name in names(contrasts_list)) {
  # Get the contrast matrix
  contrast <- contrasts_list[[contrast_name]]
  
  # Apply the contrast and perform the LRT
  lrt <- glmLRT(fit, contrast = contrast)
  
  # Extract the top genes from the LRT result
  top_genes <- topTags(lrt, n = Inf)$table
  
  # Filter significant genes based on p-value < 0.001 and log fold change (not absolute)
  significant_genes <- top_genes[top_genes$PValue < 0.001, ]
  
  # Separate upregulated and downregulated genes
  upregulated_genes <- significant_genes[significant_genes$logFC > 1, ]
  downregulated_genes <- significant_genes[significant_genes$logFC < -1, ]
  
  # Save the lists as CSV files
  upregulated_filename <- paste0(contrast_name, "_upregulated_genes.csv")
  downregulated_filename <- paste0(contrast_name, "_downregulated_genes.csv")
  
  write.csv(upregulated_genes, file = upregulated_filename, row.names = TRUE)
  write.csv(downregulated_genes, file = downregulated_filename, row.names = TRUE)
  
  # Print message indicating the lists were saved
  cat("Saved upregulated and downregulated gene lists for", contrast_name, "\n")
}


##########################

# MA Plot
# Loop through each contrast and save the MA plots
for (contrast_name in names(contrasts_list)) {
  # Get the contrast matrix
  contrast <- contrasts_list[[contrast_name]]
  
  # Apply the contrast and perform the LRT
  lrt <- glmLRT(fit, contrast = contrast)
  
  # Step 1: Create the MA plot
  plotMA(lrt, main = paste(contrast_name, "- MA Plot"), ylim = c(-5, 5))
  
  # Step 2: Save the MA plot as a PNG file
  plot_filename_ma <- paste0(contrast_name, "_MA_plot.png")
  png(plot_filename_ma)
  plotMA(lrt, main = paste(contrast_name, "- MA Plot"), ylim = c(-5, 5))
  dev.off()  # Close the PNG device
  
  # Print message indicating the plot was saved
  cat("Saved MA plot for", contrast_name, "\n")
}

# Loop through each contrast and save the MA plots
for (contrast_name in names(contrasts_list)) {
  # Get the contrast matrix
  contrast <- contrasts_list[[contrast_name]]
  
  # Apply the contrast and perform the LRT
  lrt <- glmLRT(fit, contrast = contrast)
  
  # Check if the LRT result contains the required data for plotMA
  if (!is.null(lrt$table)) {
    
    # Step 1: Create the MA plot
    plotMA(lrt, main = paste(contrast_name, "- MA Plot"), ylim = c(-5, 5))
    
    # Step 2: Save the MA plot as a PNG file
    plot_filename_ma <- paste0(contrast_name, "_MA_plot.png")
    png(plot_filename_ma)
    plotMA(lrt, main = paste(contrast_name, "- MA Plot"), ylim = c(-5, 5))
    dev.off()  # Close the PNG device
    
    # Print message indicating the plot was saved
    cat("Saved MA plot for", contrast_name, "\n")
    
  } else {
    cat("No valid data found for", contrast_name, "\n")
  }
}


library(ggplot2)

for (contrast_name in names(contrasts_list)) {
  # Get the contrast matrix
  contrast <- contrasts_list[[contrast_name]]
  
  # Apply the contrast and perform the LRT
  lrt <- glmLRT(fit, contrast = contrast)
  
  # Check if the table exists and has the necessary components
  if (!is.null(lrt$table) && all(c("logFC", "logCPM", "PValue") %in% colnames(lrt$table))) {
    
    # Extract logFC, logCPM, and PValue for plotting
    plot_data <- data.frame(logFC = lrt$table$logFC, logCPM = lrt$table$logCPM, PValue = lrt$table$PValue)
    
    # Add a column to indicate significance based on PValue threshold (0.001)
    plot_data$Significant <- ifelse(plot_data$PValue < 0.001, "Significant", "Non-significant")
    
    # Create an MA plot using ggplot2, coloring by significance
    p <- ggplot(plot_data, aes(x = logCPM, y = logFC, color = Significant)) +
      geom_point(alpha = 0.6) +
      geom_hline(yintercept = c(-1, 1), col = "red", linetype = "dashed") +
      scale_color_manual(values = c("Non-significant" = "gray", "Significant" = "blue")) +
      labs(title = paste(contrast_name, "- Significance threshold: p < 0.001"),
           x = "Average logCPM",
           y = "logFC") +
      theme_minimal() +
      ylim(-5, 5)
    
    # Print the plot
    print(p)
    
    # Save the plot
    plot_filename <- paste0(contrast_name, "_MA_plot_manual.png")
    ggsave(plot_filename, p)
    
    cat("Saved manual MA plot for", contrast_name, "\n")
    
  } else {
    cat("MA plot could not be created for", contrast_name, "due to missing or invalid data\n")
  }
}

#####

for (contrast_name in names(contrasts_list)) {
  # Get the contrast matrix
  contrast <- contrasts_list[[contrast_name]]
  
  # Apply the contrast and perform the LRT
  lrt <- glmLRT(fit, contrast = contrast)
  
  # Check if the table exists and has the necessary components
  if (!is.null(lrt$table) && all(c("logFC", "logCPM", "PValue") %in% colnames(lrt$table))) {
    
    # Extract logFC, logCPM, and PValue for plotting
    plot_data <- data.frame(logFC = lrt$table$logFC, logCPM = lrt$table$logCPM, PValue = lrt$table$PValue)
    
    # Filter genes with logFC > 1 and PValue < 0.001
    plot_data <- plot_data[abs(plot_data$logFC) > 1, ]
    
    # Add a column to indicate significance based on PValue threshold (0.001)
    plot_data$Significant <- ifelse(plot_data$PValue < 0.001, "Significant", "Non-significant")
    
    # Create an MA plot using ggplot2, coloring by significance
    p <- ggplot(plot_data, aes(x = logCPM, y = logFC, color = Significant)) +
      geom_point(alpha = 0.6) +
      geom_hline(yintercept = c(-1, 1), col = "red", linetype = "dashed") +
      scale_color_manual(values = c("Non-significant" = "gray", "Significant" = "blue")) +
      labs(title = paste(contrast_name, "- MA Plot\nSignificance Threshold: p < 0.001, |logFC| > 1"),
           x = "Average logCPM",
           y = "logFC") +
      theme_minimal() +
      ylim(-5, 5)
    
    # Print the plot
    print(p)
    
    # Save the plot
    plot_filename <- paste0(contrast_name, "_MA_plot_manual_filtered.png")
    ggsave(plot_filename, p)
    
    cat("Saved filtered MA plot for", contrast_name, "\n")
    
  } else {
    cat("MA plot could not be created for", contrast_name, "due to missing or invalid data\n")
  }
}


#####
library(pheatmap)

# Select the top 100 DEGs (or any number you want)
top_genes <- rownames(topTags(lrt, n = 100)$table)

# Assuming `metadata` contains a "Condition" column to cluster by

# Create a simplified annotation dataframe, removing the sample IDs
annotation_col <- data.frame(Condition = metadata$Condition)
rownames(annotation_col) <- rownames(metadata)  # Keep sample IDs to match the heatmap, but won't display

# Plot the heatmap
pheatmap(expression_matrix, 
         scale = "row",                         # Scale by row to Z-scores
         annotation_col = annotation_col,        # Annotate samples by condition only
         show_rownames = FALSE,                  # Hide gene names
         show_colnames = FALSE,                  # Hide sample names (IDs)
         main = "Heatmap of Top 100 DEGs Clustered by Condition")

######

# Convert the lrt$table to a data frame
lrt_table <- as.data.frame(lrt$table)

# Sort by lowest p-value and then by largest absolute log fold change (logFC)
lrt_table_sorted <- lrt_table[order(lrt_table$PValue, -abs(lrt_table$logFC)), ]

# Extract the top 25 genes based on this combined criteria
top_genes <- rownames(lrt_table_sorted)[1:25]

# Merge the top genes with gene names from `all.gene.meta`
top_genes_with_names <- merge(data.frame(ID = top_genes), all.gene.meta, by.x = "ID", by.y = "ID", all.x = TRUE)

# Check the result
head(top_genes_with_names)

# Get the gene names
top_gene_names <- top_genes_with_names$Gene.name

# Subset the expression matrix using gene IDs
top_25_expression_matrix <- dge_filtered$counts[top_genes_with_names$ID, ]

# Ensure rownames are the gene names for labeling
rownames(top_25_expression_matrix) <- top_gene_names

# Plot the heatmap with gene names
pheatmap(as.matrix(top_25_expression_matrix), 
         scale = "row",                         # Scale by row to Z-scores
         annotation_col = annotation_col,        # Annotate samples by condition
         show_rownames = TRUE,                   # Show gene names in the rows
         show_colnames = TRUE,                  # Hide sample names
         main = "Heatmap of Top Genes\nRanked by Lowest P-value and Largest LogFC")

# Extract the top 25 genes based on this combined criteria
top_genes <- rownames(lrt_table_sorted)[1:40]

# Merge the top genes with gene names from `all.gene.meta`
top_genes_with_names <- merge(data.frame(ID = top_genes), all.gene.meta, by.x = "ID", by.y = "ID", all.x = TRUE)

# Check the result
head(top_genes_with_names)

# Get the gene names
top_gene_names <- top_genes_with_names$Gene.name

# Subset the expression matrix using gene IDs
top_25_expression_matrix <- dge_filtered$counts[top_genes_with_names$ID, ]

# Ensure rownames are the gene names for labeling
rownames(top_25_expression_matrix) <- top_gene_names

# Plot the heatmap with gene names
pheatmap(as.matrix(top_25_expression_matrix), 
         scale = "row",                         # Scale by row to Z-scores
         annotation_col = annotation_col,        # Annotate samples by condition
         show_rownames = TRUE,                   # Show gene names in the rows
         show_colnames = TRUE,                  # Hide sample names
         main = "Heatmap of Top Genes\nRanked by Lowest P-value and Largest LogFC")

###################
# Boxplot of Log Counts per Million (LogCPM)

log_cpm <- cpm(dge_filtered, log = TRUE)

boxplot(log_cpm, las = 2, col = rainbow(ncol(log_cpm)),
        main = "LogCPM Distribution across Samples", ylab = "LogCPM")

library(RColorBrewer)

# Define colors by condition
condition_colors <- brewer.pal(length(unique(metadata$Condition)), "Set1")

# Match conditions to sample columns
sample_colors <- condition_colors[as.numeric(as.factor(metadata$Condition))]

boxplot(log_cpm, las = 2, col = sample_colors,
        main = "LogCPM Distribution across Samples", ylab = "LogCPM")

# Optionally, add a legend for conditions
legend("bottomleft", legend = unique(metadata$Condition), fill = condition_colors, cex = 0.8)


log_cpm <- cpm(dge_test, log = TRUE)

boxplot(log_cpm, las = 2, col = rainbow(ncol(log_cpm)),
        main = "LogCPM Distribution across Samples
        before filtering/normalisation", ylab = "LogCPM")

library(RColorBrewer)

# Define colors by condition
condition_colors <- brewer.pal(length(unique(metadata$Condition)), "Set1")

# Match conditions to sample columns
sample_colors <- condition_colors[as.numeric(as.factor(metadata$Condition))]

boxplot(log_cpm, las = 2, col = sample_colors,
        main = "LogCPM Distribution across Samples", ylab = "LogCPM")

# Optionally, add a legend for conditions
legend("bottomleft", legend = unique(metadata$Condition), fill = condition_colors, cex = 0.8)

################

# Hierarchical Clustering Dendrogram
dist_matrix <- dist(t(log_cpm))  # Distance matrix based on logCPM values
clustering <- hclust(dist_matrix)

# Plot the dendrogram
plot(clustering, labels = metadata$Condition, main = "Hierarchical Clustering Dendrogram")

#################

BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "enrichplot", "pathview"))

# Load the required packages
library(clusterProfiler)
library(org.Mm.eg.db)  # Mouse gene annotation package
library(enrichplot)    # For visualization
library(pathview)      # For KEGG pathway visualization
library(ggplot2)


###########################
# KEGG Enrichment Loop
for (contrast_name in names(contrasts_list)) {
  # Get the contrast matrix
  contrast <- contrasts_list[[contrast_name]]
  
  # Apply the contrast and perform the LRT
  lrt <- glmLRT(fit, contrast = contrast)
  
  # Extract top differentially expressed genes based on significance
  top_genes <- topTags(lrt, n = Inf)$table
  
  # Filter genes with adjusted p-value < 0.05 and log fold change > 1
  significant_genes <- top_genes[top_genes$PValue < 0.05 & abs(top_genes$logFC) > 1, ]
  
  # Check if significant genes exist
  if (nrow(significant_genes) == 0) {
    cat("No significant DEGs for contrast:", contrast_name, "\n")
    next
  }
  
  # Create a named vector of logFC for significant genes
  gene_list <- significant_genes$logFC
  names(gene_list) <- rownames(significant_genes)
  
  # Convert gene IDs to Entrez IDs for KEGG
  entrez_ids <- bitr(names(gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  
  #### KEGG Enrichment ####
  kegg_results <- enrichKEGG(
    gene = entrez_ids$ENTREZID,
    organism = "mmu",  # For mouse (Mus musculus)
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1
  )
  
  # If there are significant KEGG pathways, save the results and plot
  if (!is.null(kegg_results) && nrow(as.data.frame(kegg_results)) > 0) {
    # Save the KEGG enrichment results as a CSV file
    kegg_results_filename <- paste0(contrast_name, "_KEGG_enrichment_results.csv")
    write.csv(as.data.frame(kegg_results), file = kegg_results_filename, row.names = FALSE)
    
    # Create dot plot for KEGG enrichment results
    kegg_plot <- dotplot(kegg_results, showCategory = 20) +
      ggtitle(paste(contrast_name, "- KEGG Pathway Enrichment (p < 0.05, q < 0.1)")) +
      theme(axis.text.y = element_text(size = 5),  # Decrease the y-axis text size
            plot.title = element_text(size = 12, face = "bold")) +  # Adjust the title size
      theme(panel.grid.major.y = element_blank(),  # Removes grid lines between y-axis terms
            axis.ticks.y = element_blank(),        # Removes ticks on the y-axis
            plot.margin = margin(10, 10, 10, 20))
    
    # Save the KEGG enrichment plot
    kegg_plot_filename <- paste0(contrast_name, "_KEGG_enrichment_plot.png")
    ggsave(kegg_plot_filename, kegg_plot, height = 6, width = 8)
    
    cat("Saved KEGG plot and results for contrast:", contrast_name, "\n")
  } else {
    cat("No significant KEGG pathways found for contrast:", contrast_name, "\n")
  }
}
 
###########################

# GO Enrichment Loop
for (contrast_name in names(contrasts_list)) {
  # Get the contrast matrix
  contrast <- contrasts_list[[contrast_name]]
  
  # Apply the contrast and perform the LRT
  lrt <- glmLRT(fit, contrast = contrast)
  
  # Extract top differentially expressed genes based on significance
  top_genes <- topTags(lrt, n = Inf)$table
  
  # Filter genes with adjusted p-value < 0.05 and log fold change > 1
  significant_genes <- top_genes[top_genes$PValue < 0.05 & abs(top_genes$logFC) > 1, ]
  
  # Check if significant genes exist
  if (nrow(significant_genes) == 0) {
    cat("No significant DEGs for contrast:", contrast_name, "\n")
    next
  }
  
  # Create a named vector of logFC for significant genes
  gene_list <- significant_genes$logFC
  names(gene_list) <- rownames(significant_genes)
  
  # Convert gene IDs to Entrez IDs for GO
  entrez_ids <- bitr(names(gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  
  #### GO Enrichment ####
  go_results <- enrichGO(
    gene = entrez_ids$ENTREZID,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",  # Biological Process
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1
  )
  
  # If there are significant GO terms, save the results and plot
  if (!is.null(go_results) && nrow(as.data.frame(go_results)) > 0) {
    # Save the GO enrichment results as a CSV file
    go_results_filename <- paste0(contrast_name, "_GO_enrichment_results.csv")
    write.csv(as.data.frame(go_results), file = go_results_filename, row.names = FALSE)
    
    # Create dot plot for GO enrichment results
    go_plot <- dotplot(go_results, showCategory = 20) +
      ggtitle(paste(contrast_name, "- GO Enrichment (p < 0.05, q < 0.1)")) +
      theme(axis.text.y = element_text(size = 5),  # Decrease the y-axis text size
            plot.title = element_text(size = 12, face = "bold")) +  # Adjust the title size
      theme(panel.grid.major.y = element_blank(),  # Removes grid lines between y-axis terms
            axis.ticks.y = element_blank(),        # Removes ticks on the y-axis
            plot.margin = margin(10, 10, 10, 20))
    
    
    # Save the GO enrichment plot
    go_plot_filename <- paste0(contrast_name, "_GO_enrichment_plot.png")
    ggsave(go_plot_filename, go_plot, height = 6, width = 8)
    
    cat("Saved GO plot and results for contrast:", contrast_name, "\n")
  } else {
    cat("No significant GO terms found for contrast:", contrast_name, "\n")
  }
}
#####################3
# Required libraries
library(fgsea)
library(msigdbr)
library(tidyverse)
library(ggplot2)

# List all DEG files (already present in deg_files variable)
deg_folder <- "C:/Users/agsko/dev/pcm/EdgeR"
deg_files <- list.files(deg_folder, pattern = "DEG_20,000_genes.csv", full.names = TRUE)

# Load the C3 pathways for Mus musculus
c3_pathways <- msigdbr(species = "Mus musculus", category = "C3")

# Prepare the C3 pathways list: group by pathway names and list the Ensembl gene IDs
c3_pathways_list <- c3_pathways %>%
  dplyr::select(gs_name, ensembl_gene) %>%
  group_by(gs_name) %>%
  summarize(genes = list(ensembl_gene)) %>%
  deframe()

# Loop through each DEG file and run FGSEA analysis
for (deg_file in deg_files) {
  
  # Load the DEG data
  deg_data <- read.csv(deg_file)
  
  # Prepare the ranked list of genes using the 'ID' column (Ensembl IDs) and 'logFC' for ranking
  ranks <- deg_data %>%
    dplyr::filter(!is.na(ID), !is.na(logFC)) %>%
    dplyr::select(ID, logFC) %>%
    dplyr::distinct() %>%
    deframe()  # Convert to named vector (Ensembl IDs as names, logFC as values)
  
  # Run FGSEA using the C3 pathways
  fgsea_res <- fgsea(pathways = c3_pathways_list,
                     stats = ranks,
                     minSize = 15,  # Minimum size of pathways
                     maxSize = 500)
  
  # Order results by the adjusted p-value
  fgsea_res <- fgsea_res[order(fgsea_res$padj), ]
  
  # Remove the list columns from fgsea_res before saving
  fgsea_res_clean <- fgsea_res %>%
    dplyr::select(-leadingEdge)  # Remove the 'leadingEdge' column (which is a list)
  
  # Generate output file names based on the DEG file name
  contrast_name <- tools::file_path_sans_ext(basename(deg_file))  # Extract the contrast name from the file
  output_csv <- paste0("fgsea_results_", contrast_name, "_C3.csv")
  
  # Save the cleaned results to a CSV file
  write.csv(fgsea_res_clean, output_csv, row.names = FALSE)
  
  # Print a message indicating completion of analysis for this file
  message(paste("FGSEA analysis completed for", contrast_name, "and saved results as", output_csv, "and", output_plot))
}


######################################
# List all DEG files (already present in deg_files variable)
deg_folder <- "C:/Users/agsko/dev/pcm/EdgeR"
deg_files <- list.files(deg_folder, pattern = "DEG_20,000_genes.csv", full.names = TRUE)

# Load the C2 pathways for Mus musculus
c2_pathways <- msigdbr(species = "Mus musculus", category = "C2")

# Prepare the C2 pathways list: group by pathway names and list the Ensembl gene IDs
c2_pathways_list <- c2_pathways %>%
  dplyr::select(gs_name, ensembl_gene) %>%
  group_by(gs_name) %>%
  summarize(genes = list(ensembl_gene)) %>%
  deframe()

# Loop through each DEG file and run FGSEA analysis
for (deg_file in deg_files) {
  
  # Load the DEG data
  deg_data <- read.csv(deg_file)
  
  # Prepare the ranked list of genes using the 'ID' column (Ensembl IDs) and 'logFC' for ranking
  ranks <- deg_data %>%
    dplyr::filter(!is.na(ID), !is.na(logFC)) %>%
    dplyr::select(ID, logFC) %>%
    dplyr::distinct() %>%
    deframe()  # Convert to named vector (Ensembl IDs as names, logFC as values)
  
  # Run FGSEA using the C2 pathways
  fgsea_res <- fgsea(pathways = c2_pathways_list,
                     stats = ranks,
                     minSize = 15,  # Minimum size of pathways
                     maxSize = 500)
  
  # Order results by the adjusted p-value
  fgsea_res <- fgsea_res[order(fgsea_res$padj), ]
  
  # Remove the list columns from fgsea_res before saving
  fgsea_res_clean <- fgsea_res %>%
    dplyr::select(-leadingEdge)  # Remove the 'leadingEdge' column (which is a list)
  
  # Generate output file names based on the DEG file name
  contrast_name <- tools::file_path_sans_ext(basename(deg_file))  # Extract the contrast name from the file
  output_csv <- paste0("fgsea_results_", contrast_name, "_C2.csv")
  
  # Save the cleaned results to a CSV file
  write.csv(fgsea_res_clean, output_csv, row.names = FALSE)
  
  # Print a message indicating completion of analysis for this file
  message(paste("FGSEA analysis completed for", contrast_name, "and saved results as", output_csv, "and", output_plot))
}


###################
# List all DEG files (already present in deg_files variable)
deg_folder <- "C:/Users/agsko/dev/pcm/EdgeR"
deg_files <- list.files(deg_folder, pattern = "DEG_20,000_genes.csv", full.names = TRUE)

# Load the H pathways for Mus musculus
H_pathways <- msigdbr(species = "Mus musculus", category = "H")

# Prepare the H pathways list: group by pathway names and list the Ensembl gene IDs
H_pathways_list <- H_pathways %>%
  dplyr::select(gs_name, ensembl_gene) %>%
  group_by(gs_name) %>%
  summarize(genes = list(ensembl_gene)) %>%
  deframe()

# Loop through each DEG file and run FGSEA analysis
for (deg_file in deg_files) {
  
  # Load the DEG data
  deg_data <- read.csv(deg_file)
  
  # Prepare the ranked list of genes using the 'ID' column (Ensembl IDs) and 'logFC' for ranking
  ranks <- deg_data %>%
    dplyr::filter(!is.na(ID), !is.na(logFC)) %>%
    dplyr::select(ID, logFC) %>%
    dplyr::distinct() %>%
    deframe()  # Convert to named vector (Ensembl IDs as names, logFC as values)
  
  # Run FGSEA using the H pathways
  fgsea_res <- fgsea(pathways = H_pathways_list,
                     stats = ranks,
                     minSize = 15,  # Minimum size of pathways
                     maxSize = 500)
  
  # Order results by the adjusted p-value
  fgsea_res <- fgsea_res[order(fgsea_res$padj), ]
  
  # Remove the list columns from fgsea_res before saving
  fgsea_res_clean <- fgsea_res %>%
    dplyr::select(-leadingEdge)  # Remove the 'leadingEdge' column (which is a list)
  
  # Generate output file names based on the DEG file name
  contrast_name <- tools::file_path_sans_ext(basename(deg_file))  # Extract the contrast name from the file
  output_csv <- paste0("fgsea_results_", contrast_name, "_H.csv")
  
  # Save the cleaned results to a CSV file
  write.csv(fgsea_res_clean, output_csv, row.names = FALSE)
  
  # Print a message indicating completion of analysis for this file
  message(paste("FGSEA analysis completed for", contrast_name, "and saved results as", output_csv, "and", output_plot))
}

##########################
# List all DEG files (already present in deg_files variable)
deg_folder <- "C:/Users/agsko/dev/pcm/EdgeR"
deg_files <- list.files(deg_folder, pattern = "DEG_20,000_genes.csv", full.names = TRUE)

# Load the C5 pathways for Mus musculus
c5_pathways <- msigdbr(species = "Mus musculus", category = "C5")

# Prepare the C5 pathways list: group by pathway names and list the Ensembl gene IDs
c5_pathways_list <- c5_pathways %>%
  dplyr::select(gs_name, ensembl_gene) %>%
  group_by(gs_name) %>%
  summarize(genes = list(ensembl_gene)) %>%
  deframe()

# Loop through each DEG file and run FGSEA analysis
for (deg_file in deg_files) {
  
  # Load the DEG data
  deg_data <- read.csv(deg_file)
  
  # Prepare the ranked list of genes using the 'ID' column (Ensembl IDs) and 'logFC' for ranking
  ranks <- deg_data %>%
    dplyr::filter(!is.na(ID), !is.na(logFC)) %>%
    dplyr::select(ID, logFC) %>%
    dplyr::distinct() %>%
    deframe()  # Convert to named vector (Ensembl IDs as names, logFC as values)
  
  # Run FGSEA using the C5 pathways
  fgsea_res <- fgsea(pathways = c5_pathways_list,
                     stats = ranks,
                     minSize = 15,  # Minimum size of pathways
                     maxSize = 500)
  
  # Order results by the adjusted p-value
  fgsea_res <- fgsea_res[order(fgsea_res$padj), ]
  
  # Remove the list columns from fgsea_res before saving
  fgsea_res_clean <- fgsea_res %>%
    dplyr::select(-leadingEdge)  # Remove the 'leadingEdge' column (which is a list)
  
  # Generate output file names based on the DEG file name
  contrast_name <- tools::file_path_sans_ext(basename(deg_file))  # Extract the contrast name from the file
  output_csv <- paste0("fgsea_results_", contrast_name, "_C5.csv")
  
  # Save the cleaned results to a CSV file
  write.csv(fgsea_res_clean, output_csv, row.names = FALSE)
  
  # Print a message indicating completion of analysis for this file
  message(paste("FGSEA analysis completed for", contrast_name, "and saved results as", output_csv, "and", output_plot))
}

########################
# List all DEG files (already present in deg_files variable)
deg_folder <- "C:/Users/agsko/dev/pcm/EdgeR"
deg_files <- list.files(deg_folder, pattern = "DEG_20,000_genes.csv", full.names = TRUE)

# Load the C8 pathways for Mus musculus
c8_pathways <- msigdbr(species = "Mus musculus", category = "C8")

# Prepare the C8 pathways list: group by pathway names and list the Ensembl gene IDs
c8_pathways_list <- c8_pathways %>%
  dplyr::select(gs_name, ensembl_gene) %>%
  group_by(gs_name) %>%
  summarize(genes = list(ensembl_gene)) %>%
  deframe()

# Loop through each DEG file and run FGSEA analysis
for (deg_file in deg_files) {
  
  # Load the DEG data
  deg_data <- read.csv(deg_file)
  
  # Prepare the ranked list of genes using the 'ID' column (Ensembl IDs) and 'logFC' for ranking
  ranks <- deg_data %>%
    dplyr::filter(!is.na(ID), !is.na(logFC)) %>%
    dplyr::select(ID, logFC) %>%
    dplyr::distinct() %>%
    deframe()  # Convert to named vector (Ensembl IDs as names, logFC as values)
  
  # Run FGSEA using the C8 pathways
  fgsea_res <- fgsea(pathways = c8_pathways_list,
                     stats = ranks,
                     minSize = 15,  # Minimum size of pathways
                     maxSize = 500)
  
  # Order results by the adjusted p-value
  fgsea_res <- fgsea_res[order(fgsea_res$padj), ]
  
  # Remove the list columns from fgsea_res before saving
  fgsea_res_clean <- fgsea_res %>%
    dplyr::select(-leadingEdge)  # Remove the 'leadingEdge' column (which is a list)
  
  # Generate output file names based on the DEG file name
  contrast_name <- tools::file_path_sans_ext(basename(deg_file))  # Extract the contrast name from the file
  output_csv <- paste0("fgsea_results_", contrast_name, "_C8.csv")
  
  # Save the cleaned results to a CSV file
  write.csv(fgsea_res_clean, output_csv, row.names = FALSE)
  
  # Print a message indicating completion of analysis for this file
  message(paste("FGSEA analysis completed for", contrast_name, "and saved results as", output_csv, "and", output_plot))
}



##############


############
# Load necessary libraries
library(ggplot2)

# Load the FGSEA results for each file
fgsea_res_1 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_DMXAAvsPBS_RT_DMXAA_DEG_20,000_genes_C3.csv")
fgsea_res_2 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_DMXAAvsUre_DEG_20,000_genes_C3.csv")
fgsea_res_3 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_DMXAAvsUre_RT_DEG_20,000_genes_C3.csv")
fgsea_res_4 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RT_DMXAAvsPBS_RT_DMXAA_DEG_20,000_genes_C3.csv")
fgsea_res_5 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RT_DMXAAvsUre_DEG_20,000_genes_C3.csv")
fgsea_res_6 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RT_DMXAAvsUre_DMXAA_DEG_20,000_genes_C3.csv")
fgsea_res_7 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RT_DMXAAvsUre_RT_DEG_20,000_genes_C3.csv")
fgsea_res_8 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RTvsPBS_RT_DMXAA_DEG_20,000_genes_C3.csv")
fgsea_res_9 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RTvsUre_DEG_20,000_genes_C3.csv")
fgsea_res_10 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_UrevsPBS_RT_DMXAA_DEG_20,000_genes_C3.csv")

# Get the top 20 most enriched pathways for each fgsea result

topPathways_1 <- fgsea_res_1 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_2 <- fgsea_res_2 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_3 <- fgsea_res_3 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_4 <- fgsea_res_4 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_5 <- fgsea_res_5 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_6 <- fgsea_res_6 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_7 <- fgsea_res_7 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_8 <- fgsea_res_8 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_9 <- fgsea_res_9 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_10 <- fgsea_res_10 %>%
  dplyr::arrange(padj) %>%
  head(20)

# Plot top enriched pathways for each topPathways result

p1 <- ggplot(topPathways_1, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nregulatory gene sets (Ure_DMXAA vs PBS_RT_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),  # Adjust y-axis text size
    plot.title = element_text(size = 12)    # Adjust title text size
  )
print(p1)

p2 <- ggplot(topPathways_2, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nregulatory gene sets (Ure_DMXAA vs Ure)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p2)

p3 <- ggplot(topPathways_3, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nregulatory gene sets (Ure_DMXAA vs Ure_RT)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p3)

p4 <- ggplot(topPathways_4, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nregulatory gene sets (Ure_RT_DMXAA vs PBS_RT_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p4)

p5 <- ggplot(topPathways_5, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nregulatory gene sets (Ure_RT_DMXAA vs Ure)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p5)

p6 <- ggplot(topPathways_6, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nregulatory gene sets (Ure_RT_DMXAA vs Ure_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p6)

p7 <- ggplot(topPathways_7, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nregulatory gene sets (Ure_RT_DMXAA vs Ure_RT)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p7)

p8 <- ggplot(topPathways_8, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nregulatory gene sets (Ure_RT vs PBS_RT_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p8)

p9 <- ggplot(topPathways_9, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nregulatory gene sets (Ure_RT vs Ure)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p9)

p10 <- ggplot(topPathways_10, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nregulatory gene sets (Ure vs PBS_RT_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p10)

########################
# Load the FGSEA results for each file
fgsea_res_1 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_DMXAAvsPBS_RT_DMXAA_DEG_20,000_genes_C2.csv")
fgsea_res_2 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_DMXAAvsUre_DEG_20,000_genes_C2.csv")
fgsea_res_3 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_DMXAAvsUre_RT_DEG_20,000_genes_C2.csv")
fgsea_res_4 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RT_DMXAAvsPBS_RT_DMXAA_DEG_20,000_genes_C2.csv")
fgsea_res_5 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RT_DMXAAvsUre_DEG_20,000_genes_C2.csv")
fgsea_res_6 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RT_DMXAAvsUre_DMXAA_DEG_20,000_genes_C2.csv")
fgsea_res_7 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RT_DMXAAvsUre_RT_DEG_20,000_genes_C2.csv")
fgsea_res_8 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RTvsPBS_RT_DMXAA_DEG_20,000_genes_C2.csv")
fgsea_res_9 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RTvsUre_DEG_20,000_genes_C2.csv")
fgsea_res_10 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_UrevsPBS_RT_DMXAA_DEG_20,000_genes_C2.csv")

# Get the top 20 most enriched pathways for fgsea_res_1 to fgsea_res_10

topPathways_1 <- fgsea_res_1 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_2 <- fgsea_res_2 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_3 <- fgsea_res_3 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_4 <- fgsea_res_4 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_5 <- fgsea_res_5 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_6 <- fgsea_res_6 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_7 <- fgsea_res_7 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_8 <- fgsea_res_8 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_9 <- fgsea_res_9 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_10 <- fgsea_res_10 %>%
  dplyr::arrange(padj) %>%
  head(20)

# Plot top enriched pathways for each topPathways result

p1 <- ggplot(topPathways_1, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\ncurated gene sets (Ure_DMXAA vs PBS_RT_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 5),  # Adjust y-axis text size
    plot.title = element_text(size = 12)    # Adjust title text size
  )
print(p1)

p2 <- ggplot(topPathways_2, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\ncurated gene sets (Ure_DMXAA vs Ure)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 5),
    plot.title = element_text(size = 12)
  )
print(p2)

p3 <- ggplot(topPathways_3, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\ncurated gene sets (Ure_DMXAA vs Ure_RT)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 5),
    plot.title = element_text(size = 12)
  )
print(p3)

p4 <- ggplot(topPathways_4, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\ncurated gene sets (Ure_RT_DMXAA vs PBS_RT_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 5),
    plot.title = element_text(size = 12)
  )
print(p4)

p5 <- ggplot(topPathways_5, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\ncurated gene sets (Ure_RT_DMXAA vs Ure)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 5),
    plot.title = element_text(size = 12)
  )
print(p5)

p6 <- ggplot(topPathways_6, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\ncurated gene sets (Ure_RT_DMXAA vs Ure_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 5),
    plot.title = element_text(size = 12)
  )
print(p6)

p7 <- ggplot(topPathways_7, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\ncurated gene sets (Ure_RT_DMXAA vs Ure_RT)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 5),
    plot.title = element_text(size = 12)
  )
print(p7)

p8 <- ggplot(topPathways_8, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\ncurated gene sets (Ure_RT vs PBS_RT_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 5),
    plot.title = element_text(size = 12)
  )
print(p8)

p9 <- ggplot(topPathways_9, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\ncurated gene sets (Ure_RT vs Ure)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 5),
    plot.title = element_text(size = 12)
  )
print(p9)

p10 <- ggplot(topPathways_10, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\ncurated gene sets (Ure vs PBS_RT_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 5),
    plot.title = element_text(size = 12)
  )
print(p10)

#######################

# Load the FGSEA results for each file
fgsea_res_1 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_DMXAAvsPBS_RT_DMXAA_DEG_20,000_genes_C5.csv")
fgsea_res_2 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_DMXAAvsUre_DEG_20,000_genes_C5.csv")
fgsea_res_3 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_DMXAAvsUre_RT_DEG_20,000_genes_C5.csv")
fgsea_res_4 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RT_DMXAAvsPBS_RT_DMXAA_DEG_20,000_genes_C5.csv")
fgsea_res_5 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RT_DMXAAvsUre_DEG_20,000_genes_C5.csv")
fgsea_res_6 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RT_DMXAAvsUre_DMXAA_DEG_20,000_genes_C5.csv")
fgsea_res_7 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RT_DMXAAvsUre_RT_DEG_20,000_genes_C5.csv")
fgsea_res_8 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RTvsPBS_RT_DMXAA_DEG_20,000_genes_C5.csv")
fgsea_res_9 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RTvsUre_DEG_20,000_genes_C5.csv")
fgsea_res_10 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_UrevsPBS_RT_DMXAA_DEG_20,000_genes_C5.csv")

# Get the top 20 most enriched pathways for each fgsea result

# Get the top 20 most enriched pathways for fgsea_res_1 to fgsea_res_10

topPathways_1 <- fgsea_res_1 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_2 <- fgsea_res_2 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_3 <- fgsea_res_3 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_4 <- fgsea_res_4 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_5 <- fgsea_res_5 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_6 <- fgsea_res_6 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_7 <- fgsea_res_7 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_8 <- fgsea_res_8 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_9 <- fgsea_res_9 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_10 <- fgsea_res_10 %>%
  dplyr::arrange(padj) %>%
  head(20)

# Plot top enriched pathways for each topPathways result

p1 <- ggplot(topPathways_1, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nontology gene sets (Ure_DMXAA vs PBS_RT_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 5),  # Adjust y-axis text size
    plot.title = element_text(size = 12)    # Adjust title text size
  )
print(p1)

p2 <- ggplot(topPathways_2, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nontology gene sets (Ure_DMXAA vs Ure)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 5),
    plot.title = element_text(size = 12)
  )
print(p2)

p3 <- ggplot(topPathways_3, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nontology gene sets (Ure_DMXAA vs Ure_RT)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 5),
    plot.title = element_text(size = 12)
  )
print(p3)

p4 <- ggplot(topPathways_4, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nontology gene sets (Ure_RT_DMXAA vs PBS_RT_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 5),
    plot.title = element_text(size = 12)
  )
print(p4)

p5 <- ggplot(topPathways_5, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nontology gene sets (Ure_RT_DMXAA vs Ure)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 5),
    plot.title = element_text(size = 12)
  )
print(p5)

p6 <- ggplot(topPathways_6, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nontology gene sets (Ure_RT_DMXAA vs Ure_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 5),
    plot.title = element_text(size = 12)
  )
print(p6)

p7 <- ggplot(topPathways_7, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nontology gene sets (Ure_RT_DMXAA vs Ure_RT)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 5),
    plot.title = element_text(size = 12)
  )
print(p7)

p8 <- ggplot(topPathways_8, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nontology gene sets (Ure_RT vs PBS_RT_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 5),
    plot.title = element_text(size = 12)
  )
print(p8)

p9 <- ggplot(topPathways_9, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nontology gene sets (Ure_RT vs Ure)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 5),
    plot.title = element_text(size = 12)
  )
print(p9)

p10 <- ggplot(topPathways_10, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nontology gene sets (Ure vs PBS_RT_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 5),
    plot.title = element_text(size = 12)
  )
print(p10)

################
# Load necessary libraries
library(ggplot2)

# Load the FGSEA results for each file
fgsea_res_1 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_DMXAAvsPBS_RT_DMXAA_DEG_20,000_genes_H.csv")
fgsea_res_2 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_DMXAAvsUre_DEG_20,000_genes_H.csv")
fgsea_res_3 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_DMXAAvsUre_RT_DEG_20,000_genes_H.csv")
fgsea_res_4 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RT_DMXAAvsPBS_RT_DMXAA_DEG_20,000_genes_H.csv")
fgsea_res_5 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RT_DMXAAvsUre_DEG_20,000_genes_H.csv")
fgsea_res_6 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RT_DMXAAvsUre_DMXAA_DEG_20,000_genes_H.csv")
fgsea_res_7 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RT_DMXAAvsUre_RT_DEG_20,000_genes_H.csv")
fgsea_res_8 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RTvsPBS_RT_DMXAA_DEG_20,000_genes_H.csv")
fgsea_res_9 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RTvsUre_DEG_20,000_genes_H.csv")
fgsea_res_10 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_UrevsPBS_RT_DMXAA_DEG_20,000_genes_H.csv")

# Get the top 20 most enriched pathways for each fgsea result

topPathways_1 <- fgsea_res_1 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_2 <- fgsea_res_2 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_3 <- fgsea_res_3 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_4 <- fgsea_res_4 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_5 <- fgsea_res_5 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_6 <- fgsea_res_6 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_7 <- fgsea_res_7 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_8 <- fgsea_res_8 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_9 <- fgsea_res_9 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_10 <- fgsea_res_10 %>%
  dplyr::arrange(padj) %>%
  head(20)

# Plot top enriched pathways for each topPathways result

p1 <- ggplot(topPathways_1, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nhallmark gene sets (Ure_DMXAA vs PBS_RT_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),  # Adjust y-axis text size
    plot.title = element_text(size = 12)    # Adjust title text size
  )
print(p1)

p2 <- ggplot(topPathways_2, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nhallmark gene sets (Ure_DMXAA vs Ure)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p2)

p3 <- ggplot(topPathways_3, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nhallmark gene sets (Ure_DMXAA vs Ure_RT)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p3)

p4 <- ggplot(topPathways_4, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nhallmark gene sets (Ure_RT_DMXAA vs PBS_RT_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p4)

p5 <- ggplot(topPathways_5, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nhallmark gene sets (Ure_RT_DMXAA vs Ure)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p5)

p6 <- ggplot(topPathways_6, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nhallmark gene sets (Ure_RT_DMXAA vs Ure_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p6)

p7 <- ggplot(topPathways_7, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nhallmark gene sets (Ure_RT_DMXAA vs Ure_RT)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p7)

p8 <- ggplot(topPathways_8, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nhallmark gene sets (Ure_RT vs PBS_RT_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p8)

p9 <- ggplot(topPathways_9, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nhallmark gene sets (Ure_RT vs Ure)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p9)

p10 <- ggplot(topPathways_10, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\nhallmark gene sets (Ure vs PBS_RT_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p10)

##################
library(ggplot2)

# Load the FGSEA results for each file
fgsea_res_1 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_DMXAAvsPBS_RT_DMXAA_DEG_20,000_genes_C8.csv")
fgsea_res_2 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_DMXAAvsUre_DEG_20,000_genes_C8.csv")
fgsea_res_3 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_DMXAAvsUre_RT_DEG_20,000_genes_C8.csv")
fgsea_res_4 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RT_DMXAAvsPBS_RT_DMXAA_DEG_20,000_genes_C8.csv")
fgsea_res_5 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RT_DMXAAvsUre_DEG_20,000_genes_C8.csv")
fgsea_res_6 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RT_DMXAAvsUre_DMXAA_DEG_20,000_genes_C8.csv")
fgsea_res_7 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RT_DMXAAvsUre_RT_DEG_20,000_genes_C8.csv")
fgsea_res_8 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RTvsPBS_RT_DMXAA_DEG_20,000_genes_C8.csv")
fgsea_res_9 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_Ure_RTvsUre_DEG_20,000_genes_C8.csv")
fgsea_res_10 <- read.csv("C:/Users/agsko/dev/pcm/EdgeR/fgsea_results_UrevsPBS_RT_DMXAA_DEG_20,000_genes_C8.csv")

# Get the top 20 most enriched pathways for each fgsea result

# Get the top 20 most enriched pathways for fgsea_res_1 to fgsea_res_10

topPathways_1 <- fgsea_res_1 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_2 <- fgsea_res_2 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_3 <- fgsea_res_3 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_4 <- fgsea_res_4 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_5 <- fgsea_res_5 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_6 <- fgsea_res_6 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_7 <- fgsea_res_7 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_8 <- fgsea_res_8 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_9 <- fgsea_res_9 %>%
  dplyr::arrange(padj) %>%
  head(20)

topPathways_10 <- fgsea_res_10 %>%
  dplyr::arrange(padj) %>%
  head(20)

# Plot top enriched pathways for each topPathways result

p1 <- ggplot(topPathways_1, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\ncell type signature gene sets (Ure_DMXAA vs PBS_RT_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),  # Adjust y-axis text size
    plot.title = element_text(size = 12)    # Adjust title text size
  )
print(p1)

p2 <- ggplot(topPathways_2, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\ncell type signature gene sets (Ure_DMXAA vs Ure)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p2)

p3 <- ggplot(topPathways_3, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\ncell type signature gene sets (Ure_DMXAA vs Ure_RT)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p3)

p4 <- ggplot(topPathways_4, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\ncell type signature gene sets (Ure_RT_DMXAA vs PBS_RT_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p4)

p5 <- ggplot(topPathways_5, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\ncell type signature gene sets (Ure_RT_DMXAA vs Ure)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p5)

p6 <- ggplot(topPathways_6, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\ncell type signature gene sets (Ure_RT_DMXAA vs Ure_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p6)

p7 <- ggplot(topPathways_7, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\ncell type signature gene sets (Ure_RT_DMXAA vs Ure_RT)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p7)

p8 <- ggplot(topPathways_8, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\ncell type signature gene sets (Ure_RT vs PBS_RT_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p8)

p9 <- ggplot(topPathways_9, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\ncell type signature gene sets (Ure_RT vs Ure)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p9)

p10 <- ggplot(topPathways_10, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj)) +
  coord_flip() +
  labs(title = "Top Enriched Pathways in\ncell type signature gene sets (Ure vs PBS_RT_DMXAA)",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 12)
  )
print(p10)

#########################
#COmmon KEGG pathways

# Load necessary libraries
library(dplyr)

# Define the paths to the KEGG enrichment results files
files <- list(
  "Ure_RTvsUre" = "C:/Users/agsko/dev/pcm/EdgeR/KEGG results/Ure_RTvsUre_KEGG_enrichment_results.csv",
  "Ure_RT_DMXAAvsUre_RT" = "C:/Users/agsko/dev/pcm/EdgeR/KEGG results/Ure_RT_DMXAAvsUre_RT_KEGG_enrichment_results.csv",
  "Ure_RT_DMXAAvsUre_DMXAA" = "C:/Users/agsko/dev/pcm/EdgeR/KEGG results/Ure_RT_DMXAAvsUre_DMXAA_KEGG_enrichment_results.csv",
  "Ure_RT_DMXAAvsUre" = "C:/Users/agsko/dev/pcm/EdgeR/KEGG results/Ure_RT_DMXAA_vs_Ure_KEGG_enrichment_results.csv",
  "Ure_DMXAAvsUre" = "C:/Users/agsko/dev/pcm/EdgeR/KEGG results/Ure_DMXAAvsUre_KEGG_enrichment_results.csv"
)

# Read the files into a list of data frames
kegg_data <- lapply(files, read.csv)

# Create a column to track the source of the data
for (i in seq_along(kegg_data)) {
  kegg_data[[i]] <- kegg_data[[i]] %>% mutate(Source = names(files)[i])
}

# Combine all data into one dataframe
combined_kegg <- bind_rows(kegg_data)

# Find common occurrences based on the 'ID' column
common_kegg <- combined_kegg %>%
  group_by(ID) %>%
  filter(n() > 1) %>%
  summarise(across(everything(), ~ toString(unique(.))), .groups = "drop")

# Write the results to a new CSV file
write.csv(common_kegg, "C:/Users/agsko/dev/pcm/EdgeR/KEGG results/KEGG_results_common_occurences.csv", row.names = FALSE)

# View the resulting data
print(common_kegg)

#################

#COmmon GO terms

# Load necessary libraries
library(dplyr)

# Define the paths to the GO enrichment results files
go_files <- list(
  "Ure_RT_DMXAA_vs_Ure_RT" = "C:/Users/agsko/dev/pcm/EdgeR/GO/Ure_RT_DMXAA_vs_Ure_RT_GO_enrichment_results.csv",
  "Ure_RT_vs_Ure" = "C:/Users/agsko/dev/pcm/EdgeR/GO/Ure_RT_vs_Ure_GO_enrichment_results.csv",
  "Ure_DMXAA_vs_Ure" = "C:/Users/agsko/dev/pcm/EdgeR/GO/Ure_DMXAA_vs_Ure_GO_enrichment_results.csv",
  "Ure_RT_DMXAA_vs_Ure_DMXAA" = "C:/Users/agsko/dev/pcm/EdgeR/GO/Ure_RT_DMXAA_vs_Ure_DMXAA_GO_enrichment_results.csv",
  "Ure_RT_DMXAA_vs_Ure" = "C:/Users/agsko/dev/pcm/EdgeR/GO/Ure_RT_DMXAA_vs_Ure_GO_enrichment_results.csv"
)

# Read the files into a list of data frames
go_data <- lapply(go_files, read.csv)

# Create a column to track the source of the data
for (i in seq_along(go_data)) {
  go_data[[i]] <- go_data[[i]] %>% mutate(Source = names(go_files)[i])
}

# Combine all data into one dataframe
combined_go <- bind_rows(go_data)

# Find common occurrences based on the 'ID' column
common_go <- combined_go %>%
  group_by(ID) %>%
  filter(n() > 1) %>%
  summarise(across(everything(), ~ toString(unique(.))), .groups = "drop")

# Write the results to a new CSV file
write.csv(common_go, "C:/Users/agsko/dev/pcm/EdgeR/GO/GO_results_common_occurences.csv", row.names = FALSE)

# View the resulting data
print(common_go)

# List of GO terms to search for
go_terms <- c(
  "negative regulation of peptidase activity",
  "negative regulation of hydrolase activity",
  "immunoglobulin production",
  "adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains",
  "negative regulation of leukocyte activation",
  "humoral immune response",
  "leukocyte cell-cell adhesion",
  "regulation of leukocyte cell-cell adhesion"
)

# File paths
file_paths <- list(
  "Ure_RT_DMXAA_vs_Ure_RT" = "C:/Users/agsko/dev/pcm/EdgeR/GO/Ure_RT_DMXAA_vs_Ure_RT_GO_enrichment_results.csv",
  "Ure_RT_vs_Ure" = "C:/Users/agsko/dev/pcm/EdgeR/GO/Ure_RT_vs_Ure_GO_enrichment_results.csv",
  "Ure_DMXAA_vs_Ure" = "C:/Users/agsko/dev/pcm/EdgeR/GO/Ure_DMXAA_vs_Ure_GO_enrichment_results.csv",
  "Ure_RT_DMXAA_vs_Ure_DMXAA" = "C:/Users/agsko/dev/pcm/EdgeR/GO/Ure_RT_DMXAA_vs_Ure_DMXAA_GO_enrichment_results.csv",
  "Ure_RT_DMXAA_vs_Ure" = "C:/Users/agsko/dev/pcm/EdgeR/GO/Ure_RT_DMXAA_vs_Ure_GO_enrichment_results.csv"
)

# Initialize an empty data frame to collect all the results
final_results <- data.frame()

# Function to filter the file based on the GO terms
filter_go_terms <- function(file_path) {
  # Read the file
  df <- read.csv(file_path)
  
  # Filter the rows that contain the specified GO terms
  filtered_df <- df[df$Description %in% go_terms, ]
  
  # Select the required columns
  filtered_df <- filtered_df[, c("Count", "Description", "FoldEnrichment", "p.adjust")]
  
  return(filtered_df)
}

# Loop through all files and filter based on GO terms
for (file_name in names(file_paths)) {
  file_path <- file_paths[[file_name]]
  
  # Apply the filtering function
  filtered_data <- filter_go_terms(file_path)
  
  # If any data found, add it to the final results with a column specifying the file it came from
  if (nrow(filtered_data) > 0) {
    filtered_data$Source <- file_name
    final_results <- rbind(final_results, filtered_data)
  }
}

# Save the final results to a new file
write.csv(final_results, "C:/Users/agsko/dev/pcm/EdgeR/GO/filtered_GO_terms.csv", row.names = FALSE)



"Ure_RTvsUre" = "C:/Users/agsko/dev/pcm/EdgeR/KEGG results/Ure_RTvsUre_KEGG_enrichment_results.csv",
"Ure_RT_DMXAAvsUre_RT" = "C:/Users/agsko/dev/pcm/EdgeR/KEGG results/Ure_RT_DMXAAvsUre_RT_KEGG_enrichment_results.csv",
"Ure_RT_DMXAAvsUre_DMXAA" = "C:/Users/agsko/dev/pcm/EdgeR/KEGG results/Ure_RT_DMXAAvsUre_DMXAA_KEGG_enrichment_results.csv",
"Ure_RT_DMXAAvsUre" = "C:/Users/agsko/dev/pcm/EdgeR/KEGG results/Ure_RT_DMXAA_vs_Ure_KEGG_enrichment_results.csv",
"Ure_DMXAAvsUre" = "C:/Users/agsko/dev/pcm/EdgeR/KEGG results/Ure_DMXAAvsUre_KEGG_enrichment_results.csv"