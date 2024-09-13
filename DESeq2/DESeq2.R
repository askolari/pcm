# If DESeq2 is not installed, you can install it with:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

library(DESeq2)
library(readxl)


# Load the raw counts data
raw_counts_file <- "C:/Users/agsko/dev/pcm/raw_counts.csv"
raw_counts <- read.csv(raw_counts_file)

# Load the colData
coldata_file <- "C:/Users/agsko/dev/pcm/coldata.xlsx"
coldata <- read_excel(coldata_file)

# Rename the first column of raw counts to "Gene_ID" for clarity
colnames(raw_counts)[1] <- "Gene_ID"

# Store the Gene_ID and Gene.name columns as a new data frame
gene_info <- data.frame(
  Gene_ID = raw_counts$Gene_ID,
  Gene_Name = raw_counts$Gene.name
)

# Check the first few rows to ensure the mapping is correct
head(gene_info)

# Match the sample names between raw_counts and coldata
sample_columns <- colnames(raw_counts)[colnames(raw_counts) %in% coldata$Sample]

# Ensure that the sample names in raw_counts match the Sample column in coldata
if (!all(sample_columns %in% coldata$Sample)) {
  stop("Sample names in raw_counts do not match coldata!")
}

# Check the column names to ensure everything is correct
colnames(raw_counts)

# Confirm dimensions match
if (length(sample_columns) != nrow(coldata)) {
  stop("Number of samples in colData does not match number of sample columns in raw counts!")
}

# Subset raw_counts to include only Gene_ID and the sample columns
count_matrix <- as.matrix(raw_counts[, sample_columns])

# Ensure the rownames of the count matrix are the gene IDs
rownames(count_matrix) <- raw_counts$Gene_ID

# Convert the 'Condition' column in coldata to a factor
coldata$Condition <- as.factor(coldata$Condition)

# Check structure of coldata and ensure conditions match what you expect
str(coldata)


# Create the DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~ Condition)

# Run the DESeq function
dds <- DESeq(dds)

# Load the necessary libraries
library(DESeq2)

# Ensure the DESeq2 object (dds) has already been created and DESeq has been run
dds <- DESeq(dds)

# Define all contrasts to analyze
contrast_list <- list(
  Ure_vs_Ure_DMXAA = c("Condition", "Ure", "Ure_DMXAA"),
  Ure_vs_Ure_RT = c("Condition", "Ure", "Ure_RT"),
  Ure_vs_Ure_RT_DMXAA = c("Condition", "Ure", "Ure_RT_DMXAA"),
  Ure_vs_PBS_RT_DMXAA = c("Condition", "Ure", "PBS_RT_DMXAA"),
  Ure_DMXAA_vs_Ure_RT = c("Condition", "Ure_DMXAA", "Ure_RT"),
  Ure_DMXAA_vs_Ure_RT_DMXAA = c("Condition", "Ure_DMXAA", "Ure_RT_DMXAA"),
  Ure_DMXAA_vs_PBS_RT_DMXAA = c("Condition", "Ure_DMXAA", "PBS_RT_DMXAA"),
  Ure_RT_vs_Ure_RT_DMXAA = c("Condition", "Ure_RT", "Ure_RT_DMXAA"),
  Ure_RT_vs_PBS_RT_DMXAA = c("Condition", "Ure_RT", "PBS_RT_DMXAA"),
  Ure_RT_DMXAA_vs_PBS_RT_DMXAA = c("Condition", "Ure_RT_DMXAA", "PBS_RT_DMXAA")
)

# Loop over each contrast
for (contrast_name in names(contrast_list)) {
  # Get the contrast
  contrast <- contrast_list[[contrast_name]]
  
  # Run DESeq results for this contrast
  results_contrast <- results(dds, contrast = contrast)
  
  # Filter out rows with NA in padj or log2FoldChange
  filtered_results <- results_contrast[!is.na(results_contrast$padj) & !is.na(results_contrast$log2FoldChange), ]
  
  # Now apply the filtering criteria for significance
  significant_genes <- filtered_results[filtered_results$padj < 0.05 & abs(filtered_results$log2FoldChange) > 1, ]
  
  # Print the number of significant genes for the contrast
  cat("Number of significant genes for", contrast_name, ":", nrow(significant_genes), "\n")
  
  # Specify the file paths where you want to save the significant genes
  output_csv <- paste0("significant_genes_", contrast_name, ".csv")
  
  # Save the significant genes to CSV and Excel files
  write.csv(significant_genes, file = output_csv, row.names = TRUE)
  
  # Save the full results for the contrast (including non-significant genes)
  full_csv <- paste0("results_", contrast_name, ".csv")
  write.csv(as.data.frame(results_contrast), full_csv, row.names = TRUE)
}

#Number of significant genes for Ure_vs_Ure_DMXAA : 85 
#Number of significant genes for Ure_vs_Ure_RT : 219 
#Number of significant genes for Ure_vs_Ure_RT_DMXAA : 857 
#Number of significant genes for Ure_vs_PBS_RT_DMXAA : 445 
#Number of significant genes for Ure_DMXAA_vs_Ure_RT : 181 
#Number of significant genes for Ure_DMXAA_vs_Ure_RT_DMXAA : 534 
#Number of significant genes for Ure_DMXAA_vs_PBS_RT_DMXAA : 530 
#Number of significant genes for Ure_RT_vs_Ure_RT_DMXAA : 19 
#Number of significant genes for Ure_RT_vs_PBS_RT_DMXAA : 0 
#Number of significant genes for Ure_RT_DMXAA_vs_PBS_RT_DMXAA : 24


# Log2 transform the baseMean values
log_baseMean_values <- log2(baseMean_values + 1)  # Adding 1 to avoid log(0)
hist(log_baseMean_values, breaks = 50, main = "Distribution of log2(baseMean) across all conditions", 
     xlab = "log2(BaseMean)", ylab = "Frequency")

# Check the metadata of the dds object
mcols(dds)

setwd("C:/Users/agsko/dev/pcm/DESeq2")
getwd()

install.packages("VennDiagram")
install.packages("readr")

# Load necessary libraries
library(VennDiagram)
library(readr)

library(VennDiagram)

# Load necessary libraries
library(readr)
library(dplyr)

# List of file paths for each contrast
file_paths <- list(
  Ure_vs_Ure_DMXAA = "C:/Users/agsko/dev/pcm/DESeq2/significant_genes_Ure_vs_Ure_DMXAA.csv",
  Ure_vs_Ure_RT = "C:/Users/agsko/dev/pcm/DESeq2/significant_genes_Ure_vs_Ure_RT.csv",
  Ure_vs_Ure_RT_DMXAA = "C:/Users/agsko/dev/pcm/DESeq2/significant_genes_Ure_vs_Ure_RT_DMXAA.csv",
  Ure_DMXAA_vs_PBS_RT_DMXAA = "C:/Users/agsko/dev/pcm/DESeq2/significant_genes_Ure_DMXAA_vs_PBS_RT_DMXAA.csv",
  Ure_DMXAA_vs_Ure_RT = "C:/Users/agsko/dev/pcm/DESeq2/significant_genes_Ure_DMXAA_vs_Ure_RT.csv",
  Ure_DMXAA_vs_Ure_RT_DMXAA = "C:/Users/agsko/dev/pcm/DESeq2/significant_genes_Ure_DMXAA_vs_Ure_RT_DMXAA.csv",
  Ure_RT_DMXAA_vs_PBS_RT_DMXAA = "C:/Users/agsko/dev/pcm/DESeq2/significant_genes_Ure_RT_DMXAA_vs_PBS_RT_DMXAA.csv",
  Ure_RT_vs_PBS_RT_DMXAA = "C:/Users/agsko/dev/pcm/DESeq2/significant_genes_Ure_RT_vs_PBS_RT_DMXAA.csv",
  Ure_RT_vs_Ure_RT_DMXAA = "C:/Users/agsko/dev/pcm/DESeq2/significant_genes_Ure_RT_vs_Ure_RT_DMXAA.csv",
  Ure_vs_PBS_RT_DMXAA = "C:/Users/agsko/dev/pcm/DESeq2/significant_genes_Ure_vs_PBS_RT_DMXAA.csv"
)

# Initialize an empty list to store processed data
contrast_data <- list()

# Load all files and handle unnamed first column
for (contrast in names(file_paths)) {
  # Read the data, assume first column is unnamed
  data <- read_csv(file_paths[[contrast]], col_names = TRUE)
  
  # Rename the first unnamed column to 'GeneID'
  colnames(data)[1] <- "GeneID"
  
  # Print the first few rows of the data to confirm
  cat("First few rows for", contrast, ":\n")
  print(head(data))
  
  # Store the processed data in the list
  contrast_data[[contrast]] <- data
}

# The contrast_data list now holds the cleaned data for each contrast.
# For example, access the first contrast with: contrast_data$Ure_vs_Ure_DMXAA



