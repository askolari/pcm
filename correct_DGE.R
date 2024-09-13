
# Preliminary Setup

install.packages("knitr")

# Load required packages
library(knitr)
library(here)

knitr::opts_chunk$set(echo = TRUE, retina = 2, cache = TRUE, autodep = TRUE,
                      fig.height = 8, fig.width = 8,
                      cache.path = paste0(
                        here::here("cache", "rnaseq-cluster"),
                        .Platform$file.sep))

# Install CRAN packages
install.packages(c("pheatmap", "here", "ggplot2", "dplyr", "tibble", "glue"))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("edgeR", "limma", "DESeq2", "matrixStats"))
ainstall.packages("stringr")
install.packages("openxlsx")

setwd("C:/Users/agsko/dev/pcm")

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

setwd("C:/Users/agsko/dev/pcm")


set.seed(1986)

# Data Loading and Preprocessing

## Loading Data

```{r load_rna_data}
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


# Load necessary libraries
library(SummarizedExperiment)
library(openxlsx)


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

# Load necessary libraries
install.packages("dplyr")
library(dplyr)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("S4Vectors")

library(S4Vectors)

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

BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)


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
dge <- calcNormFactors(dge)

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

############################################

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

# Perform voom transformation
v <- voom(dge_filtered, design, plot = TRUE)

# Extract the voom-transformed expression values (logCPM)
voom_logcpm <- v$E

# Fit the linear model
fit <- lmFit(voom_logcpm, design)

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


