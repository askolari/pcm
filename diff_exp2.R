# Define a new contrast comparing RT-treated conditions with non-RT conditions
contrast_rt_vs_non_rt <- makeContrasts(
  RT_vs_NonRT = (Ure_RT + Ure_RT_DMXAA + PBS_RT_DMXAA) / 3 - (Ure + Ure_DMXAA) / 2,
  levels = design
)

# Fit the linear model and apply contrasts
fit_rt_vs_non_rt <- contrasts.fit(fit, contrast_rt_vs_non_rt)
fit_rt_vs_non_rt <- eBayes(fit_rt_vs_non_rt)

# Get the top differentially expressed genes for RT vs Non-RT
top_genes_rt_vs_non_rt <- topTable(fit_rt_vs_non_rt, coef = "RT_vs_NonRT", number = Inf, p.value = 0.05, lfc = 1)

# Save the results to a CSV file
write.csv(top_genes_rt_vs_non_rt, file = "RT_vs_NonRT_DEG_results.csv", row.names = TRUE)

# Visualize the top genes using a volcano plot
library(ggplot2)
top_genes_rt_vs_non_rt$Significant <- ifelse(top_genes_rt_vs_non_rt$adj.P.Val < 0.05 & abs(top_genes_rt_vs_non_rt$logFC) > 1, "Yes", "No")

volcano_plot <- ggplot(top_genes_rt_vs_non_rt, aes(x = logFC, y = -log10(P.Value), color = Significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot: RT vs Non-RT",
       x = "Log Fold Change",
       y = "-log10(P-Value)") +
  theme_minimal()

# Print the volcano plot
print(volcano_plot)

# Create a heatmap for the differentially expressed genes
rt_vs_non_rt_genes <- rownames(top_genes_rt_vs_non_rt)

# Extract top differentially expressed genes for RT vs Non-RT
top_genes_rt_vs_non_rt <- topTable(fit_rt_vs_non_rt, coef = "RT_vs_NonRT", number = Inf, p.value = 0.05, lfc = 1)

# Check if any genes were found
if (nrow(top_genes_rt_vs_non_rt) == 0) {
  stop("No differentially expressed genes found.")
}


# Functional enrichment for RT vs Non-RT genes
ego_rt_vs_non_rt <- enrichGO(gene = rownames(top_genes_rt_vs_non_rt), 
                             OrgDb = org.Mm.eg.db, 
                             keyType = "ENSEMBL",
                             ont = "BP",  # Biological Process
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.1)

dotplot(ego_rt_vs_non_rt, showCategory = 10) +
  ggtitle("GO Enrichment for RT vs Non-RT Genes") +
  theme(axis.text.y = element_text(size = 5))

# KEGG pathway enrichment analysis
# Perform KEGG pathway enrichment analysis
kegg_enrich_rt_vs_non_rt <- enrichKEGG(gene = entrez_ids$ENTREZID, 
                                       organism = 'mmu',  # 'mmu' for mouse
                                       pvalueCutoff = 0.05)

# Check if any pathways were enriched
if (!is.null(kegg_enrich_rt_vs_non_rt) && nrow(kegg_enrich_rt_vs_non_rt) > 0) {
  # Visualize the top KEGG pathways
  dotplot(kegg_enrich_rt_vs_non_rt, showCategory = 10) +
    ggtitle("KEGG Pathway Enrichment for RT vs Non-RT Genes") +
    theme(axis.text.y = element_text(size = 5))
} else {
  cat("No significant KEGG pathways found.\n")
}

library(fgsea)

# Ensure you have the ranked gene list, e.g., ranked_gene_list <- sort(...)


# Load necessary libraries
library(fgsea)
library(msigdbr)
library(ggplot2)

# Load MSigDB Hallmark gene sets (replace with your msigdb_gene_sets if already loaded)
msigdb_gene_sets <- msigdbr(species = "Mus musculus", category = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Ensure you have the ranked gene list
# Example: ranked_gene_list <- sort(your_gene_list, decreasing = TRUE)
# Make sure to have a named numeric vector where names are gene symbols or IDs and values are ranking scores.

# Load necessary library
library(dplyr)

# Step 1: Create a ranking score for each gene
# Combine logFC and p-value into a ranking score
# Rank = logFC * -log10(P.Value)
dge_results <- read.csv("Ure_vs_PBS_RT_DMXAA_DEG_results.csv", row.names = 1)
significant_genes <- significant_genes %>%
  mutate(RankingScore = logFC * -log10(P.Value))

# Step 2: Create a named vector for GSEA input
# The names should be gene identifiers, and the values should be the ranking scores
ranked_gene_list <- significant_genes$RankingScore
names(ranked_gene_list) <- rownames(significant_genes)

# Step 3: Sort the ranked gene list in decreasing order
ranked_gene_list <- sort(ranked_gene_list, decreasing = TRUE)

# Preview the top of the ranked gene list
head(ranked_gene_list)

# Perform GSEA using the Hallmark gene sets
fgsea_res <- fgsea(pathways = msigdb_gene_sets, stats = ranked_gene_list)

# Filter for significant results (adjusted p-value < 0.05)
fgsea_res <- fgsea_res[fgsea_res$padj < 0.05, ]

# View the top significant pathways
head(fgsea_res)
#Empty data.table (0 rows and 8 cols): pathway,pval,padj,log2err,ES,NES...

# Plot enrichment for a specific pathway (e.g., the first significant pathway)
plotEnrichment(pathway = msigdb_gene_sets[[fgsea_res$pathway[1]]], stats = ranked_gene_list) +
  ggtitle(paste("Enrichment Plot for", fgsea_res$pathway[1]))


# Plot enrichment for a specific pathway
plotEnrichment(pathway = msigdb_gene_sets$h[[1]], stats = ranked_gene_list)

# Dot plot for significant pathways
plotGseaTable(fgsea_res, pathways = names(fgsea_res), stats = ranked_gene_list)

