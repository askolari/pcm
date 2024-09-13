# Set appropriate column names for the design matrix
colnames(design) <- levels(metadata_combined$Condition)

# Check the design matrix
head(design)

# Define contrasts for the comparisons you are interested in
contrasts <- makeContrasts(
  Ure_vs_Ure_RT = Ure - Ure_RT,
  Ure_vs_Ure_DMXAA = Ure - Ure_DMXAA,
  Ure_DMXAA_vs_Ure_RT_DMXAA = Ure_DMXAA - Ure_RT_DMXAA,
  Ure_vs_PBS_RT_DMXAA = Ure - PBS_RT_DMXAA,
  Ure_RT_DMXAA_vs_PBS_RT_DMXAA = Ure_RT_DMXAA - PBS_RT_DMXAA,
  Ure_vs_PBS_RT_DMXAA = Ure - PBS_RT_DMXAA,
  levels = design
)

# Fit the linear model
fit <- lmFit(voom_logcpm, design)

# Apply the contrasts to the fitted model
fit_contrasts <- contrasts.fit(fit, contrasts)

# Perform empirical Bayes moderation
fit_contrasts <- eBayes(fit_contrasts)

# Check the top differentially expressed genes for each contrast
topTable(fit_contrasts, coef = "Ure_vs_Ure_RT")

# Define a list of contrast names
contrast_names <- c("Ure_vs_Ure_RT", "Ure_vs_Ure_DMXAA", 
                    "Ure_DMXAA_vs_Ure_RT_DMXAA", "Ure_vs_PBS_RT_DMXAA", 
                    "Ure_RT_DMXAA_vs_PBS_RT_DMXAA")

# Loop through each contrast and save the results
for (contrast in contrast_names) {
  # Extract top differentially expressed genes for the contrast
  top_genes <- topTable(fit_contrasts, coef = contrast, number = Inf)
  
  # Save to a CSV file
  output_filename <- paste0(contrast, "_DEG_results.csv")
  write.csv(top_genes, file = output_filename, row.names = TRUE)
  
  # Print message to indicate saving progress
  cat("Saved results for", contrast, "to", output_filename, "\n")
}

library(ggplot2)


for (contrast in contrast_names) {
  # Extract top differentially expressed genes for the contrast
  top_genes <- topTable(fit_contrasts, coef = contrast, number = Inf)
  
  # Check if top_genes is not empty
  if (nrow(top_genes) == 0) {
    cat("No genes found for contrast:", contrast, "\n")
    next  # Skip to the next iteration if no genes are found
  }
  
  # Add a column for significance
  top_genes$Significant <- ifelse(top_genes$adj.P.Val < 0.05 & abs(top_genes$logFC) > 1, "Yes", "No")
  
  # Debug: Print out a summary of the data
  cat("Contrast:", contrast, "\n")
  print(head(top_genes))
  
  # Create the volcano plot
  p <- ggplot(top_genes, aes(x = logFC, y = -log10(P.Value), color = Significant)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("grey", "red")) +
    labs(title = paste0("Volcano Plot: ", contrast),
         x = "Log Fold Change",
         y = "-log10(P-Value)") +
    theme_minimal()
  
  # Print the plot
  print(p)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")

#Gene Ontology Enrichment Analysis
#Starting with the whole dataset

library(clusterProfiler)
library(org.Mm.eg.db) 

gene_list <- top_genes$logFC
names(gene_list) <- rownames(top_genes)

ego <- enrichGO(gene = names(gene_list[gene_list > 1]), 
                OrgDb = org.Mm.eg.db, 
                keyType = "ENSEMBL",
                ont = "BP",  # Biological Process
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.1)

dotplot(ego, showCategory = 10)

library(ggplot2)

# Adjust the label size using theme
dotplot(ego) + 
  theme(axis.text.y = element_text(size = 5))  # You can change 8 to a smaller or larger value

#Specific for each contrast
# Extract top genes from a specific contrast, e.g., "Ure_vs_Ure_RT"
top_genes <- topTable(fit_contrasts, coef = "Ure_vs_Ure_RT", number = Inf)
top_genes <- top_genes[top_genes$adj.P.Val < 0.05 & abs(top_genes$logFC) > 1, ]

gene_list <- top_genes$logFC
names(gene_list) <- rownames(top_genes)

ego <- enrichGO(gene = names(gene_list), 
                OrgDb = org.Mm.eg.db, 
                keyType = "ENSEMBL",
                ont = "BP",  # Biological Process
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.1)
dotplot(ego, showCategory = 10)

head(as.data.frame(ego))
{[1] ID          Description GeneRatio   BgRatio    
  [5] pvalue      p.adjust    qvalue      geneID     
  [9] Count      
  <0 rows> (or 0-length row.names)}

ego <- enrichGO(gene = names(gene_list),
                OrgDb = org.Mm.eg.db,
                keyType = "ENSEMBL",
                ont = "BP",  # Biological Process
                pvalueCutoff = 0.14,  # Increase the cutoff to see more terms
                qvalueCutoff = 0.1)

dotplot(ego, showCategory = 10) +
  theme(axis.text.y = element_text(size = 5))

library(clusterProfiler)
library(org.Mm.eg.db)  # Assuming you're working with mouse data

# List of contrast names
contrast_names <- c("Ure_vs_Ure_RT", "Ure_vs_Ure_DMXAA", 
                    "Ure_DMXAA_vs_Ure_RT_DMXAA", "Ure_vs_PBS_RT_DMXAA", 
                    "Ure_RT_DMXAA_vs_PBS_RT_DMXAA")

# Loop through each contrast
for (contrast in contrast_names) {
  
  # Extract top differentially expressed genes for the contrast
  top_genes <- topTable(fit_contrasts, coef = contrast, number = Inf)
  
  # Filter significant genes based on adjusted p-value and log fold change
  top_genes <- top_genes[top_genes$adj.P.Val < 0.05 & abs(top_genes$logFC) > 1, ]
  
  # Create a named vector for log fold changes
  gene_list <- top_genes$logFC
  names(gene_list) <- rownames(top_genes)
  
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
      ggtitle(paste0("GO Enrichment for ", contrast)) +
      theme(axis.text.y = element_text(size = 5))
    
    # Display the plot
    print(plot)
  } else {
    cat("No significant GO terms found for contrast:", contrast, "\n")
  }
}

#GSEA Analysis
# Initialize an empty list to store gene rankings
combined_gene_list <- list()

for (contrast in contrast_names) {
  # Extract the top genes for the contrast
  top_genes <- topTable(fit_contrasts, coef = contrast, number = Inf)
  
  # Append the logFC values and gene names to the combined list
  combined_gene_list[[contrast]] <- setNames(top_genes$logFC, rownames(top_genes))
}

# Combine all contrasts into one ranked list
combined_gene_list <- unlist(combined_gene_list)

# Rank the combined gene list
combined_gene_list <- sort(combined_gene_list, decreasing = TRUE)

# Ensure that you have installed the fgsea and msigdbr packages
library(fgsea)
library(msigdbr)

# Load MSigDB pathways for mouse (C2 category)
pathways <- msigdbr(species = "Mus musculus", category = "C2") %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Run GSEA using fgseaMultilevel (recommended)
fgsea_res <- fgseaMultilevel(pathways = pathways, stats = combined_gene_list)

# Filter significant results
fgsea_res <- fgsea_res[fgsea_res$padj < 0.05, ]

# View top significant results
head(fgsea_res)

# Plot top pathway enrichment (replace 'pathways[[1]]' with a specific pathway if needed)
plotEnrichment(pathway = pathways[[1]], stats = combined_gene_list)

# Plot the distribution of log fold changes in the combined gene list
hist(combined_gene_list, breaks = 50, main = "Distribution of Log Fold Changes", 
     xlab = "Log Fold Change", col = "skyblue", border = "white")

# Initialize the combined gene list and names list
combined_gene_list <- c()
combined_gene_names <- c()

for (contrast in contrast_names) {
  top_genes <- topTable(fit_contrasts, coef = contrast, number = Inf)
  
  # Create a new ranking score based on logFC and p-value
  combined_score <- top_genes$logFC * -log10(top_genes$P.Value)
  
  # Append the scores and gene names to the combined lists
  combined_gene_list <- c(combined_gene_list, combined_score)
  combined_gene_names <- c(combined_gene_names, rownames(top_genes))
}

# Assign names to the combined gene list after the loop to ensure lengths match
names(combined_gene_list) <- combined_gene_names

# Inspect the top genes by the new ranking criteria
head(sort(combined_gene_list, decreasing = TRUE))

# Re-run fgsea with the new combined gene list
fgsea_res <- fgseaMultilevel(
  pathways = pathways, 
  stats = combined_gene_list, 
  minSize = 10,    # Adjust the pathway size thresholds as needed
  maxSize = 500
)

# Filter for significant results
fgsea_res <- fgsea_res[fgsea_res$padj < 0.05, ]

# View the top significant results
head(fgsea_res)

# Remove duplicate gene names by keeping the highest score
unique_combined_gene_list <- tapply(combined_gene_list, names(combined_gene_list), max)

# Re-run fgsea with the deduplicated gene list
fgsea_res <- fgseaMultilevel(
  pathways = pathways, 
  stats = unique_combined_gene_list, 
  minSize = 100,    # Adjust the pathway size thresholds as needed
  maxSize = 10000
)

# Filter for significant results
fgsea_res <- fgsea_res[fgsea_res$padj < 0.05, ]

# View the top significant results
head(fgsea_res)
#Empty data.table (0 rows and 8 cols): pathway,pval,padj,log2err,ES,NES...

########################################################

# Initialize a table to store the categorization results for each gene
categorization_table <- data.frame(GeneID = rownames(voom_logcpm), 
                                   First_Change_Condition = NA, 
                                   First_Change_Direction = NA, 
                                   Category = "Static", stringsAsFactors = FALSE)

# Define the baseline condition (PBS_RT_DMXAA)
baseline_contrast <- "Ure_vs_PBS_RT_DMXAA"

# Define the contrasts to compare against the baseline
contrast_names <- c("Ure_vs_Ure_RT", "Ure_vs_Ure_DMXAA", 
                    "Ure_DMXAA_vs_Ure_RT_DMXAA", "Ure_vs_PBS_RT_DMXAA", 
                    "Ure_RT_DMXAA_vs_PBS_RT_DMXAA")

# Loop through each contrast and update the categorization table
for (contrast in contrast_names) {
  # Extract top genes for the contrast
  top_genes <- topTable(fit_contrasts, coef = contrast, number = Inf)
  
  # Identify significant genes based on adjusted p-value and log fold change thresholds
  significant_genes <- top_genes[top_genes$adj.P.Val < 0.05 & abs(top_genes$logFC) > 1, ]
  
  # Update the categorization table for significant genes
  for (gene in rownames(significant_genes)) {
    if (is.na(categorization_table[categorization_table$GeneID == gene, "First_Change_Condition"])) {
      # Update the first significant change condition and direction
      categorization_table[categorization_table$GeneID == gene, "First_Change_Condition"] <- contrast
      categorization_table[categorization_table$GeneID == gene, "First_Change_Direction"] <- 
        ifelse(significant_genes[gene, "logFC"] > 0, "Upregulated", "Downregulated")
      
      # Assign categories based on response
      categorization_table[categorization_table$GeneID == gene, "Category"] <- ifelse(
        contrast == baseline_contrast, "Baseline Response", "Condition-Specific Response"
      )
    }
  }
}

# Print the categorization results
head(categorization_table)

# Count the number of genes in each category
table(categorization_table$Category)

library(ggplot2)

# Visualize the distribution of genes across categories
ggplot(categorization_table, aes(x = Category)) +
  geom_bar(fill = "steelblue") +
  labs(title = "Gene Categorization Based on Response to Conditions",
       x = "Category",
       y = "Number of Genes") +
  theme_minimal()

# Create a heatmap for genes with Baseline Response
baseline_genes <- categorization_table[categorization_table$Category == "Baseline Response", "GeneID"]

# Subset the logCPM data for the baseline response genes
baseline_logCPM <- voom_logcpm[baseline_genes, ]

# Plot a heatmap of the baseline response genes
pheatmap(baseline_logCPM,
         annotation_col = metadata_combined[c("Condition", "RT_status")],
         show_rownames = FALSE,
         show_colnames = TRUE,
         cluster_cols = TRUE, 
         clustering_distance_rows = "correlation",
         color = heat_pal)

library(clusterProfiler)
library(org.Mm.eg.db)

# Subset genes for Baseline Response and Condition-Specific Response
baseline_genes <- categorization_table[categorization_table$Category == "Baseline Response", "GeneID"]
condition_specific_genes <- categorization_table[categorization_table$Category == "Condition-Specific Response", "GeneID"]

# Functional enrichment for Baseline Response genes
ego_baseline <- enrichGO(gene = baseline_genes, 
                         OrgDb = org.Mm.eg.db, 
                         keyType = "ENSEMBL",
                         ont = "BP",  # Biological Process
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.1)

# Functional enrichment for Condition-Specific Response genes
ego_condition_specific <- enrichGO(gene = condition_specific_genes, 
                                   OrgDb = org.Mm.eg.db, 
                                   keyType = "ENSEMBL",
                                   ont = "BP",  # Biological Process
                                   pvalueCutoff = 0.05,
                                   qvalueCutoff = 0.1)

# Visualize the results using dot plots
# Baseline Response Enrichment
dotplot(ego_baseline, showCategory = 10) +
  ggtitle("GO Enrichment for Baseline Response Genes") +
  theme(axis.text.y = element_text(size = 5))

# Condition-Specific Response Enrichment
dotplot(ego_condition_specific, showCategory = 10) +
  ggtitle("GO Enrichment for Condition-Specific Response Genes") +
  theme(axis.text.y = element_text(size = 5))

# If you want to save the enrichment results
baseline_results <- as.data.frame(ego_baseline)
write.csv(baseline_results, "Baseline_Response_GO_Enrichment.csv", row.names = FALSE)

condition_specific_results <- as.data.frame(ego_condition_specific)
write.csv(condition_specific_results, "Condition_Specific_Response_GO_Enrichment.csv", row.names = FALSE)


# Extract condition-specific genes
condition_specific_genes <- categorization_table[categorization_table$Category == "Condition-Specific Response", "GeneID"]

# Subset the expression data for these genes
condition_specific_expression <- voom_logcpm[condition_specific_genes, ]

# Calculate the mean expression of each gene in each condition
condition_means <- rowMeans(condition_specific_expression)

# Create a data frame to store the gene and its associated condition
gene_condition_mapping <- data.frame(GeneID = condition_specific_genes, AssociatedCondition = NA)

# Loop through each gene and find the condition with the highest mean expression
for (gene in condition_specific_genes) {
  # Find the condition with the highest mean expression for this gene
  condition_means <- colMeans(condition_specific_expression[gene, , drop = FALSE])
  associated_condition <- names(which.max(condition_means))
  
  # Store the associated condition
  gene_condition_mapping[gene_condition_mapping$GeneID == gene, "AssociatedCondition"] <- associated_condition
}

# View the gene-to-condition mapping
head(gene_condition_mapping)

library(pheatmap)

# Visualize the expression patterns of condition-specific response genes
pheatmap(condition_specific_expression,
         annotation_col = metadata_combined[c("Condition", "RT_status")],
         show_rownames = FALSE,
         show_colnames = TRUE,
         cluster_cols = TRUE,
         clustering_distance_rows = "correlation",
         color = heat_pal,
         breaks = seq(-6, 6, length.out = length(heat_pal) + 1),
         main = "Expression Patterns of Condition-Specific Response Genes")

# Correlation analysis between gene expression and conditions
correlation_results <- data.frame(GeneID = condition_specific_genes, Correlation = NA)

for (gene in condition_specific_genes) {
  # Extract expression values for the gene
  expression_values <- condition_specific_expression[gene, ]
  
  # Perform correlation analysis (e.g., Pearson correlation)
  correlation <- cor(expression_values, as.numeric(factor(metadata_combined$Condition)), method = "pearson")
  
  # Store the correlation result
  correlation_results[correlation_results$GeneID == gene, "Correlation"] <- correlation
}

# View the correlation results
head(correlation_results)

view(correlation_results)
# Save the correlation results to a CSV file
write.csv(correlation_results, file = "condition_specific_gene_correlations.csv", row.names = FALSE)

install.packages("Hmisc")
library(Hmisc)  # For calculating correlations with p-values

# Initialize a data frame to store correlation results and p-values
correlation_results_with_pvals <- data.frame(GeneID = condition_specific_genes, Correlation = NA, P.Value = NA)
# Apply the Benjamini-Hochberg (BH) adjustment for multiple testing correction
correlation_results_with_pvals$Adj.P.Value <- p.adjust(correlation_results_with_pvals$P.Value, method = "BH")

# Check the first few rows to confirm the adjustment worked
head(correlation_results_with_pvals)


# Loop through each gene and calculate the Pearson correlation and p-value
for (gene in condition_specific_genes) {
  # Extract expression values for the gene
  expression_values <- condition_specific_expression[gene, ]
  
  # Calculate the Pearson correlation and p-value
  cor_test <- rcorr(expression_values, as.numeric(factor(metadata_combined$Condition)), type = "pearson")
  
  # Extract the correlation coefficient and p-value
  correlation <- cor_test$r[1, 2]
  p_value <- cor_test$P[1, 2]
  
  # Store the results
  correlation_results_with_pvals[correlation_results_with_pvals$GeneID == gene, "Correlation"] <- correlation
  correlation_results_with_pvals[correlation_results_with_pvals$GeneID == gene, "P.Value"] <- p_value
}

# Apply the Benjamini-Hochberg (BH) correction for multiple testing
correlation_results_with_pvals$Adj.P.Value <- p.adjust(correlation_results_with_pvals$P.Value, method = "BH")

# View the results
head(correlation_results_with_pvals)

# Save the correlation results with p-values and adjusted p-values to a CSV file
write.csv(correlation_results_with_pvals, file = "condition_specific_gene_correlations_with_pvals.csv", row.names = FALSE)

# Filter for significant correlations based on adjusted p-value (e.g., FDR < 0.05)
significant_correlations <- correlation_results_with_pvals[correlation_results_with_pvals$Adj.P.Value < 0.5, ]

# View the significant correlations
head(significant_correlations)

# Save the significant correlations to a CSV file
write.csv(significant_correlations, file = "significant_condition_specific_gene_correlations.csv", row.names = FALSE)

hist(correlation_results_with_pvals$P.Value, breaks = 50, main = "Distribution of Raw P-Values", col = "skyblue")
correlation_results_with_pvals$Adj.P.Value <- as.numeric(correlation_results_with_pvals$Adj.P.Value)

hist(correlation_results_with_pvals$Adj.P.Value, breaks = 50, main = "Distribution of Adjusted P-Values", col = "lightgreen")

###################################
#KEGG pathway analysis

library(clusterProfiler)
library(org.Mm.eg.db)

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

# Save the KEGG results to a CSV file
kegg_results <- as.data.frame(kegg_enrich)
write.csv(kegg_results, file = "kegg_enrichment_results.csv", row.names = FALSE)

cat("KEGG pathway enrichment results saved to kegg_enrichment_results.csv\n")

####################
#ReactomePA

# Install ReactomePA if not already installed
if (!requireNamespace("ReactomePA", quietly = TRUE)) {
  BiocManager::install("ReactomePA")
}

library(ReactomePA)

# Load necessary libraries
library(ReactomePA)
library(ggplot2)

# Perform Reactome pathway enrichment analysis
reactome_enrich <- enrichPathway(gene = entrez_ids$ENTREZID, 
                                 organism = "mouse", 
                                 pvalueCutoff = 0.05)

# Visualize the top enriched pathways (e.g., top 20)
dotplot(reactome_enrich, showCategory = 15) +
  ggtitle("Top Enriched Reactome Pathways") +
  theme(axis.text.y = element_text(size = 4))

# Convert the enrichment results to a data frame
reactome_results_df <- as.data.frame(reactome_enrich)

# View the first few rows of the data
head(reactome_results_df)

# Save the results to a CSV file for further exploration
write.csv(reactome_results_df, "Reactome_Pathway_Enrichment_Results.csv", row.names = FALSE)


