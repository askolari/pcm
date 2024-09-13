install.packages("WGCNA") 
library(WGCNA)
BiocManager::install("impute")
library(impute)
BiocManager::install("tidyverse")
library(tidyverse)
install.packages("magrittr")
library(magrittr)
library(tidyr)
setwd("C:/Users/agsko/dev/pcm/WGCNA")


data <- readr::read_csv("raw_counts.csv")

# ==== Load necessary libraries ====

library(DESeq2) 

# ==== Load and clean data ====
data <- readr::read_csv("C:/Users/agsko/dev/pcm/raw_counts.csv")  

# Inspect the first few rows and structure of the data
data[1:5, 1:10]        # Look at first 5 rows and 10 columns
str(data)              # Look at the structure of the dataset

# Rename the first column to "GeneId" (assuming it's gene names or IDs)
names(data)[1] <- "GeneId"
names(data)            # Look at the column names

# Transpose data: Since WGCNA expects genes as columns and samples as rows, we transpose it.
data_t <- data %>%
  tibble::column_to_rownames("GeneId") %>%  # Set 'GeneId' as rownames before transposing
  t() %>%
  as.data.frame()

# Convert the transposed data back to a tibble for better handling
data_t <- data_t %>%
  rownames_to_column(var = "Sample")

# Clean the dataset and convert to long format for visualization
col_sel <- names(data_t)[-1]  # Get all but first column name (sample names)

mdata <- data_t %>%
  tidyr::pivot_longer(
    cols = all_of(col_sel),
    names_to = "GeneId",        # Gene IDs will now be in a column called 'GeneId'
    values_to = "Expression"    # Expression values go into 'Expression' column
  ) %>%
  mutate(
    group = gsub("-.*","", Sample) %>% gsub("[.].*","", .)   # Extract treatment groups from sample names
  )


# Preview the cleaned and tidied data
head(mdata)

# Create a named vector to map Sample.IDs to groups
sample_to_group <- c(
  "X88" = "Ure", "X89" = "Ure", "X90" = "Ure",
  "X92" = "Ure_DMXAA", "X93" = "Ure_DMXAA", "X94" = "Ure_DMXAA",
  "X97" = "Ure_RT", "X99" = "Ure_RT", "X101" = "Ure_RT",
  "X103" = "Ure_RT_DMXAA", "X104" = "Ure_RT_DMXAA", "X105" = "Ure_RT_DMXAA",
  "X108" = "PBS_RT_DMXAA", "X109" = "PBS_RT_DMXAA", "X110" = "PBS_RT_DMXAA"
)

# Apply this mapping to the 'group' column
mdata <- mdata %>%
  mutate(group = sample_to_group[Sample])

# Preview the updated data
head(mdata)


# Plot groups to identify outliers
p <- mdata %>%
  ggplot(aes(x = Sample, y = as.numeric(Expression))) +   # x = Sample, y = RNA Seq count
  geom_violin() +                                         # violin plot to show distribution
  geom_point(alpha = 0.2) +                               # scatter plot
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)                # Rotate sample text
  ) +
  labs(x = "Sample Groups", y = "RNA Seq Counts") +
  facet_grid(cols = vars(group), drop = TRUE, scales = "free_x")  # Facet by groups

print(p)


# ==== Normalize counts with DESeq2 ====

de_input <- as.matrix(data[, -c(1, ncol(data))])  # Remove the first (GeneId) and last (Gene.name) columns

# Set 'GeneId' as row names
row.names(de_input) <- data$GeneId

head(de_input)

# Create metadata for samples
meta_df <- data.frame(Sample = colnames(de_input)) %>%
  mutate(Group = sample_to_group[Sample])  # Map the Sample column to Group using the sample_to_group vector
print(meta_df)

# Create DESeq2 dataset and normalize

dds <- DESeqDataSetFromMatrix(round(de_input),
                              colData=meta_df,
                              design = ~Group)

dds <- DESeq(dds)


vsd <- varianceStabilizingTransformation(dds)

library(genefilter)

wpn_vsd <- getVarianceStabilizedData(dds)
rv_wpn <- rowVars(wpn_vsd)
summary(rv_wpn)

#> summary(rv_wpn)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000000 0.000000 0.003752 0.011244 0.008197 3.287844 

q75_wpn <- quantile( rowVars(wpn_vsd), .75)  # <= original
q95_wpn <- quantile( rowVars(wpn_vsd), .95)  # <= changed to 95 quantile to reduce dataset
expr_normalized <- wpn_vsd[ rv_wpn > q95_wpn, ]

expr_normalized[1:5,1:10]

expr_normalized_df <- data.frame(expr_normalized) %>%
  mutate(
    Gene_id = row.names(expr_normalized)
  ) %>%
  pivot_longer(-Gene_id)

expr_normalized_df %>% ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    title = "Normalized and 95 quantile Expression",
    x = "treatment",
    y = "normalized expression"
  )

#Network construction
input_mat = t(expr_normalized)
input_mat[1:5,1:10]

#Group data in a dendogram to check outliers
sampleTree = hclust(dist(input_mat), method = "average")

par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Remove X94 from the dataset
input_mat_no_outliers <- input_mat[-which(rownames(input_mat) == "X94"), ]

# Re-run the WGCNA workflow on the filtered data

###############################

library(WGCNA)
allowWGCNAThreads(nThreads = 4)     

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  blockSize = 30,
  powerVector = powers,
  verbose = 5
)

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red", pos = 3, offset = 0.1
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

##############Choose power, based on observation 6,7,8

picked_power = 8
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

cor <- temp_cor     # Return cor function to original namespace

# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

netwk$colors[netwk$blockGenes[[1]]]
table(netwk$colors)


#0    1    2    3    4    5    6 
#260 1192  530  200  146   76   62 

#Relate modules (clusters) assignments to treatment groups
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

module_df[1:5,]
library(readr)

write_delim(module_df,
            file = "gene_modules.txt",
            delim = "\t")

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")

#### So brown strongly correlated with X94 alone and blue, grey and turquoise are interesting
# pick out a few modules of interest here
modules_of_interest = c("blue", "green", "brown", "red", "yellow", "turquoise")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id
# Get normalized expression for those genes
expr_normalized[1:5,1:10]

subexpr = expr_normalized[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

module_colors <- c(
  "blue" = "blue",
  "green"= "green",
  "brown" = "brown",
  "red"="red",
  "yellow"="yellow",
  "turquoise" = "turquoise"
)

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalised expression") +
  scale_color_manual(values = module_colors)


#############Generate and Export Networks
genes_of_interest = module_df %>%
  subset(colors %in% modules_of_interest)

expr_of_interest = expr_normalized[genes_of_interest$gene_id,]

expr_of_interest[1:5,1:5]

TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = picked_power)

row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)

edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )

head(edge_list)

# Export Network file to be read into Cytoscape, VisANT, etc
write_delim(edge_list,
            file = "edgelist_all.tsv",
            delim = "\t")


############################### Subsetting network by weight
library(ggplot2)

# Load the edgelist file
network <- read.table("edgelist_all.tsv", header = TRUE)

# Check the structure of the edgelist
str(network)

# Add a weight column based on the absolute value of the correlation
network$weight <- abs(network$correlation)

# Plot a histogram of the absolute correlation (weight)
ggplot(network, aes(x = weight)) +
  geom_histogram(bins = 30, fill = "blue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribution of Edge Weights (Absolute Correlations)",
       x = "Absolute Correlation (Weight)",
       y = "Frequency")


# Alternatively, you can plot a density plot
ggplot(network, aes(x = weight)) +
  geom_density(fill = "blue", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Density Plot of Edge Weights (Absolute Correlations)",
       x = "Absolute Correlation (Weight)",
       y = "Density")

#####very skewed distribution

# Subset network by keeping only edges with correlation (weights) above 0.7
edge_list_subset <- edge_list %>%
  filter(correlation > 0.7)

edge_list_subset_2 <- edge_list %>%
  filter(correlation > 0.5)

write_delim(edge_list_subset,
            file = "edgelist_subset.tsv",
            delim = "\t")

# Save the edge list for Cytoscape (gene-gene interaction list)
write.csv(edge_list_subset_2, "cytoscape_edge_list_corr_0.5.csv", row.names = FALSE)
write.csv(edge_list_subset, "cytoscape_edge_list_corr_0.7.tsv", row.names = FALSE)

# Create the igraph object from the edge list
library(igraph)

# Convert the subsetted edge list into an igraph object
g <- graph_from_data_frame(edge_list_subset, directed = FALSE)

# Calculate centrality measures
degree_centrality <- degree(g)
betweenness_centrality <- betweenness(g)
closeness_centrality <- closeness(g)

# Combine centrality values into a data frame
centrality_df <- data.frame(
  gene = names(degree_centrality),
  degree = degree_centrality,
  betweenness = betweenness_centrality,
  closeness = closeness_centrality
)

# Sort by degree to get top hub genes
centrality_df <- centrality_df[order(-centrality_df$degree), ]
head(centrality_df)  # View top hub genes based on degree centrality

# Assuming you already have module colors assigned
# Calculate module eigengenes
MEList <- moduleEigengenes(input_mat, colors = mergedColors)
MEs <- MEList$eigengenes

# Calculate the correlation matrix of module eigengenes
ME_corr <- cor(MEs)

# Plot the eigengene network heatmap
labeledHeatmap(
  Matrix = ME_corr,
  xLabels = names(MEs),
  yLabels = names(MEs),
  colorLabels = TRUE,
  colors = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Eigengene Network Heatmap"
)


# Calculate module membership (kME) for each gene
module_membership <- signedKME(input_mat, MEs)

# Combine module membership data with gene IDs
membership_df <- data.frame(
  gene = rownames(module_membership),
  module_membership
)

# Combine membership data with centrality info (if needed)
node_list <- membership_df %>%
  left_join(centrality_df, by = "gene")  # Make sure centrality_df has the 'gene' column

# Save the node list for Cytoscape
write.csv(node_list, "cytoscape_node_list_corr_0.7.csv", row.names = FALSE)


# Sort by module membership in each module (adjust the module of interest)
membership_turq_df <- membership_df[order(-membership_df$kMEturquoise), ]  # Example: for 'turquoise' module
head(membership_turq_df)  # View genes with the highest membership in a module

membership_blue_df <- membership_df[order(-membership_df$kMEblue), ]  # Example: for 'turquoise' module
head(membership_blue_df)

membership_grey_df <- membership_df[order(-membership_df$kMEgrey), ]  # Example: for 'turquoise' module
head(membership_grey_df)

membership_brown_df <- membership_df[order(-membership_df$kMEbrown), ]  # Example: for 'turquoise' module
head(membership_brown_df)

# Assuming you have calculated centrality measures (e.g., degree centrality) in a dataframe called centrality_df

# Filter membership_df to only include genes in the edge_list_subset
filtered_membership_df <- membership_df %>%
  filter(gene %in% edge_list_subset$gene1 | gene %in% edge_list_subset$gene2)

# Check the filtered data
nrow(filtered_membership_df)  # Should match the number of genes in centrality_df

# Filter based on module membership (kME) only
membership_turq_df <- membership_df %>%
  filter(kMEturquoise > 0.8)  # Use your chosen module and kME threshold

# View the top genes with the highest module membership
head(membership_turq_df)

# Save the filtered membership data (turquoise module in this case)
write.csv(membership_turq_df, "membership_turquoise_hub_genes.csv", row.names = FALSE)

# Filter based on module membership (kME) only
membership_green_df <- membership_df %>%
  filter(kMEgreen > 0.8)  # Use your chosen module and kME threshold

# View the top genes with the highest module membership
head(membership_green_df)

# Save the filtered membership data (turquoise module in this case)
write.csv(membership_green_df, "membership_green_hub_genes.csv", row.names = FALSE)

# Filter based on module membership (kME) only
membership_red_df <- membership_df %>%
  filter(kMEred > 0.8)  # Use your chosen module and kME threshold

# View the top genes with the highest module membership
head(membership_red_df)

# Save the filtered membership data (turquoise module in this case)
write.csv(membership_red_df, "membership_red_hub_genes.csv", row.names = FALSE)

# Filter based on module membership (kME) only
membership_yellow_df <- membership_df %>%
  filter(kMEyellow > 0.8)  # Use your chosen module and kME threshold

# View the top genes with the highest module membership
head(membership_yellow_df)

# Save the filtered membership data (turquoise module in this case)
write.csv(membership_yellow_df, "membership_yellow_hub_genes.csv", row.names = FALSE)

# Filter based on module membership (kME) only
membership_blue_df <- membership_df %>%
  filter(kMEblue > 0.8)  # Use your chosen module and kME threshold

# Save the filtered membership data (turquoise module in this case)
write.csv(membership_blue_df, "membership_blue_hub_genes.csv", row.names = FALSE)

# View the top genes with the highest module membership
head(membership_blue_df)

# Filter based on module membership (kME) only
membership_grey_df <- membership_df %>%
  filter(kMEgrey > 0.8)  # Use your chosen module and kME threshold

# Save the filtered membership data (turquoise module in this case)
write.csv(membership_grey_df, "membership_grey_hub_genes.csv", row.names = FALSE)

# View the top genes with the highest module membership
head(membership_grey_df)

# Filter based on module membership (kME) only
membership_brown_df <- membership_df %>%
  filter(kMEbrown > 0.8)  # Use your chosen module and kME threshold

# Save the filtered membership data (turquoise module in this case)
write.csv(membership_brown_df, "membership_brown_hub_genes.csv", row.names = FALSE)

# Filter based on module membership (kME) only
membership_brown_df <- membership_df %>%
  filter(kMEbrown > 0.8)  # Use your chosen module and kME threshold

# Save the filtered membership data (turquoise module in this case)
write.csv(membership_brown_df, "membership_brown_hub_genes.csv", row.names = FALSE)


# View the top genes with the highest module membership
head(membership_brown_df)

#############
# Filter edge list for brown module interactions
brown_edges <- edge_list %>%
  filter(module1 == "brown" & module2 == "brown")

# Save edge list for Cytoscape
write.csv(brown_edges, "cytoscape_brown_module_edge_list.csv", row.names = FALSE)

# Filter the edge list for strong interactions (correlation > 0.7 or 0.8)
brown_edge_strong <- brown_edges %>%
  filter(correlation > 0.7)  # or change to 0.8 for stronger connections

write.csv(brown_edge_strong, "cytoscape_brown_module_edge_list_filtered.csv", row.names = FALSE)

#################
# Example gene significance calculation for each module
# Create a dataframe where each sample is linked to its group (condition)
trait_data <- data.frame(
  Sample = names(sample_to_group),
  Group = sample_to_group
)

# Convert the treatment groups into a binary matrix (1 for presence, 0 for absence)
trait_matrix <- model.matrix(~ 0 + Group, data = trait_data)
colnames(trait_matrix) <- gsub("Group", "", colnames(trait_matrix))

# Make sure the order of samples matches with input_mat
rownames(trait_matrix) <- trait_data$Sample
trait_matrix <- trait_matrix[rownames(input_mat), ]

# Recalculate the module eigengenes
MEList <- moduleEigengenes(input_mat, colors = mergedColors)
MEs <- MEList$eigengenes

# Correlate module eigengenes with the traits
module_trait_cor <- cor(MEs, trait_matrix, use = "p")
module_trait_pvalues <- corPvalueStudent(module_trait_cor, nrow(input_mat))

# Plot the heatmap of module-trait relationships
labeledHeatmap(Matrix = module_trait_cor,
               xLabels = colnames(trait_matrix),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = signif(module_trait_cor, 2),
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = "Module-trait relationships")

############## URe condition
# Calculate gene significance for each module (replace with your trait of interest)
gene_significance <- cor(input_mat, trait_matrix[, "Ure"])  # Example using the 'Ure' group
module_colors <- c("blue", "brown", "turquoise", "grey")

# Get gene significance for modules of interest
gene_significance_df <- data.frame(
  module = module_df$colors,
  gene_significance = gene_significance
) %>%
  filter(module %in% module_colors)

ggplot(gene_significance_df, aes(x = module, y = gene_significance, fill = module)) +
  stat_summary(fun = mean, geom = "bar", position = "dodge") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", position = "dodge", width = 0.25) +
  scale_fill_manual(values = c("blue", "brown", "turquoise", "grey")) +
  labs(title = "Gene Significance across Modules", y = "Gene Significance", x = "Module") +
  theme_minimal()

# module_membership contains kME for all genes in each module

# Subset the membership and gene significance for the turquoise module
blue_df <- data.frame(
  gene_id = rownames(membership_df),
  gene_significance = gene_significance,
  kMEturquoise = membership_df$kMEturquoise
) %>%
  filter(module_df$colors == "blue")

# Plot module membership vs gene significance
ggplot(blue_df, aes(x = kMEturquoise, y = gene_significance)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(title = "Ure Condition: Module Membership vs. Gene Significance (Turquoise Module)",
       x = "Module Membership (kME)", y = "Gene Significance") +
  theme_minimal()

# Extract the sample names for "Ure"
ure_samples <- names(sample_to_group[sample_to_group == "Ure"])

# Subset input_mat for the "Ure" condition (rows)
ure_expr <- input_mat[ure_samples, ]

# Check the dimensions and the first few rows of ure_expr
dim(ure_expr)
head(ure_expr)

# Create a binary vector for the Ure condition
ure_group <- ifelse(sample_to_group == "Ure", 1, 0)

# Calculate gene significance as the correlation between each gene's expression and the Ure group
gene_significance_ure <- apply(input_mat, 2, function(gene_expr) cor(gene_expr, ure_group))

# Convert to a data frame for easier handling
gs_ure_df <- data.frame(
  gene_id = colnames(input_mat),
  gene_significance = gene_significance_ure,
  condition = "Ure"
)

# Create a binary vector for Ure_RT_DMXAA condition
ure_rt_dmxaa_group <- ifelse(sample_to_group == "Ure_RT_DMXAA", 1, 0)

# Calculate gene significance for Ure_RT_DMXAA
gene_significance_ure_rt_dmxaa <- apply(input_mat, 2, function(gene_expr) cor(gene_expr, ure_rt_dmxaa_group))

# Convert to data frame
gs_ure_rt_dmxaa_df <- data.frame(
  gene_id = colnames(input_mat),
  gene_significance = gene_significance_ure_rt_dmxaa,
  condition = "Ure_RT_DMXAA"
)

ure_rt_group <- ifelse(sample_to_group == "Ure_RT", 1, 0)

gene_significance_ure_rt <- apply(input_mat, 2, function(gene_expr) cor(gene_expr, ure_rt_group))

gs_ure_rt_df <- data.frame(
  gene_id = colnames(input_mat),
  gene_significance = gene_significance_ure_rt,
  condition = "Ure_RT"
)

ure_dmxaa_group <- ifelse(sample_to_group == "Ure_DMXAA", 1, 0)

gene_significance_ure_dmxaa <- apply(input_mat, 2, function(gene_expr) cor(gene_expr, ure_dmxaa_group))

gs_ure_dmxaa_df <- data.frame(
  gene_id = colnames(input_mat),
  gene_significance = gene_significance_ure_dmxaa,
  condition = "Ure_DMXAA"
)

pbs_rt_dmxaa_group <- ifelse(sample_to_group == "PBS_RT_DMXAA", 1, 0)

gene_significance_pbs_rt_dmxaa <- apply(input_mat, 2, function(gene_expr) cor(gene_expr, pbs_rt_dmxaa_group))

gs_pbs_rt_dmxaa_df <- data.frame(
  gene_id = colnames(input_mat),
  gene_significance = gene_significance_pbs_rt_dmxaa,
  condition = "PBS_RT_DMXAA"
)

# Combine all the results into one data frame
combined_gene_significance_df <- rbind(
  gs_ure_df,
  gs_ure_rt_dmxaa_df,
  gs_ure_rt_df,
  gs_ure_dmxaa_df,
  gs_pbs_rt_dmxaa_df
)

# Preview the combined results
head(combined_gene_significance_df)

combined_gene_significance_df$module <- module_df[combined_gene_significance_df$gene_id, "colors"]

# Plot gene significance across modules and conditions
ggplot(combined_gene_significance_df, aes(x = module, y = gene_significance, fill = condition)) +
  geom_boxplot() +
  facet_wrap(~condition, scales = "free_y") +
  labs(title = "Gene Significance Across Modules for Different Conditions", 
       x = "Module", 
       y = "Gene Significance") +
  theme_minimal()

PBS_RT_DMXAA_group <- ifelse(sample_to_group == "PBS_RT_DMXAA", 1, 0)

gene_significance_pbs_rt_dmxaa <- apply(input_mat, 2, function(gene_expr) cor(gene_expr, PBS_RT_DMXAA_group))


# 1. Calculate p-values for each condition and module
p_value_df <- combined_gene_significance_df %>%
  group_by(condition, module) %>%
  summarise(
    p_value = t.test(gene_significance)$p.value  # Using t-test to check if the mean significance is different from 0
  ) %>%
  ungroup()

# 2. Create a string for display with module and p-value
p_value_df <- p_value_df %>%
  mutate(
    p_label = paste0("p = ", format(p_value, digits = 3))
  )

# 3. Combine the p-value text with the original gene significance data
plot_data <- combined_gene_significance_df %>%
  left_join(p_value_df, by = c("condition", "module"))

# 4. Plot gene significance with p-value subtitles
ggplot(plot_data, aes(x = module, y = gene_significance, fill = condition)) +
  geom_boxplot() +
  facet_wrap(~condition, scales = "free") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.text = element_text(size = 10, face = "bold")  # Adjust facet labels size
  ) +
  labs(
    title = "Gene Significance Across Modules for Different Conditions",
    y = "Gene Significance",
    x = "Module"
  ) +
  scale_fill_manual(values = c(
    "Ure" = "red", 
    "Ure_RT" = "cyan", 
    "Ure_DMXAA" = "green", 
    "Ure_RT_DMXAA" = "purple", 
    "PBS_RT_DMXAA" = "orange"
  )) +
  geom_text(
    data = p_value_df,
    aes(x = module, y = Inf, label = p_label),  # Placing p-values at the top of each plot
    inherit.aes = FALSE,
    vjust = 1.5,  # Adjust the vertical placement
    size = 1,
    color = "black"
  )

############
# List of significant modules
significant_modules <- c("blue", "green", "brown", "red", "yellow", "turquoise")

# Subset the edge list by these modules
subset_edge_list <- edge_list %>%
  filter(module1 %in% significant_modules | module2 %in% significant_modules)

# Preview the subsetted edge list
head(subset_edge_list)

# Extract genes from the subsetted edge list
genes_in_modules <- unique(c(subset_edge_list$gene1, subset_edge_list$gene2))

# Separate genes by module for further analysis
genes_by_module <- subset_edge_list %>%
  filter(module1 %in% significant_modules | module2 %in% significant_modules) %>%
  distinct(gene1, module1) %>%
  arrange(module1)

# Preview genes in modules
head(genes_by_module)

write.csv(genes_by_module, file = "genes_by_module_0.7.csv", row.names = FALSE)

library(clusterProfiler)
library(org.Mm.eg.db)

# Convert Ensembl IDs (gene1) to Entrez IDs
genes_by_module <- genes_by_module %>%
  mutate(EntrezID = mapIds(org.Mm.eg.db, 
                           keys = gene1, 
                           column = "ENTREZID", 
                           keytype = "ENSEMBL", 
                           multiVals = "first"))

# Remove rows with missing Entrez IDs
genes_by_module <- genes_by_module %>%
  filter(!is.na(EntrezID))

# Now perform GO enrichment with the Entrez IDs
go_results <- enrichGO(gene = genes_by_module$EntrezID,
                       OrgDb = org.Mm.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",  # Can be "BP", "MF", "CC"
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)


# View the updated data
head(genes_by_module)
# View the results
head(go_results)

library(enrichplot)

# Filter the genes for the "blue" module
blue_genes_with_entrez <- genes_by_module %>%
  filter(module1 == "blue" & !is.na(EntrezID))

# Extract the Entrez IDs from the filtered genes (in this case, blue module)
entrez_ids <- unlist(blue_genes_with_entrez$EntrezID)

# Perform GO enrichment analysis
go_results <- enrichGO(
  gene         = entrez_ids,            # Entrez IDs of genes of interest
  OrgDb        = org.Mm.eg.db,          # Mouse gene annotation database
  keyType      = "ENTREZID",            # Key type for gene IDs
  ont          = "BP",                  # Ontology: Biological Process (BP), Molecular Function (MF), or Cellular Component (CC)
  pAdjustMethod = "BH",                 # Method for p-value adjustment
  pvalueCutoff = 0.05,                  # p-value cutoff
  qvalueCutoff = 0.2,                   # q-value cutoff
  readable     = TRUE                   # Converts Entrez IDs to gene symbols
)

barplot(go_results, showCategory = 15, title = "Top GO Terms for Blue Module")+
  theme(axis.text.y = element_text(size = 8))

dotplot(go_results, showCategory = 15, title = "Dotplot for Blue Module GO Terms")+
  theme(axis.text.y = element_text(size = 8))

# Save gene set with Entrez IDs and gene symbols for Cytoscape
# Flatten the EntrezID list into a character string for the blue module
blue_genes_with_entrez <- blue_genes_with_entrez %>%
  mutate(EntrezID = sapply(EntrezID, function(x) paste(x, collapse = ",")))

# Save gene set with Entrez IDs and gene symbols for Cytoscape
write.csv(blue_genes_with_entrez, file = "blue_genes_with_entrez.csv", row.names = FALSE)

# Filter the genes for the "turquoise" module
turquoise_genes_with_entrez <- genes_by_module %>%
  filter(module1 == "turquoise" & !is.na(EntrezID))

# Extract the Entrez IDs from the filtered genes (in this case, turquoise module)
entrez_ids_turquoise <- unlist(turquoise_genes_with_entrez$EntrezID)

# Perform GO enrichment analysis
go_results_turquoise <- enrichGO(
  gene         = entrez_ids_turquoise,
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable     = TRUE
)

barplot(go_results_turquoise, showCategory = 15, title = "Top GO Terms for Turquoise Module")+
  theme(axis.text.y = element_text(size = 6))

dotplot(go_results_turquoise, showCategory = 15, title = "Dotplot for Turquoise Module GO Terms")+
  theme(axis.text.y = element_text(size = 8))

# Save gene set with Entrez IDs and gene symbols for Cytoscape
# Flatten the EntrezID list into a character string
turquoise_genes_with_entrez <- turquoise_genes_with_entrez %>%
  mutate(EntrezID = sapply(EntrezID, function(x) paste(x, collapse = ",")))

# Save gene set with Entrez IDs and gene symbols for Cytoscape
write.csv(turquoise_genes_with_entrez, file = "turquoise_genes_with_entrez.csv", row.names = FALSE)

# Filter the genes for the "brown" module
brown_genes_with_entrez <- genes_by_module %>%
  filter(module1 == "brown" & !is.na(EntrezID))

# Extract the Entrez IDs from the filtered genes (in this case, brown module)
entrez_ids_brown <- unlist(brown_genes_with_entrez$EntrezID)

# Perform GO enrichment analysis
go_results_brown <- enrichGO(
  gene         = entrez_ids_brown,
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable     = TRUE
)

barplot(go_results_brown, showCategory = 15, title = "Top GO Terms for Brown Module")+
  theme(axis.text.y = element_text(size = 8))

dotplot(go_results_brown, showCategory = 15, title = "Dotplot for Brown Module GO Terms")+
  theme(axis.text.y = element_text(size = 8))

# Save gene set with Entrez IDs and gene symbols for Cytoscape
# Flatten the EntrezID list into a character string for the brown module
brown_genes_with_entrez <- brown_genes_with_entrez %>%
  mutate(EntrezID = sapply(EntrezID, function(x) paste(x, collapse = ",")))

# Save gene set with Entrez IDs and gene symbols for Cytoscape
write.csv(brown_genes_with_entrez, file = "brown_genes_with_entrez.csv", row.names = FALSE)

# Filter the genes for the "green" module
green_genes_with_entrez <- genes_by_module %>%
  filter(module1 == "green" & !is.na(EntrezID))

# Extract the Entrez IDs from the filtered genes (in this case, green module)
entrez_ids_green <- unlist(green_genes_with_entrez$EntrezID)

# Perform GO enrichment analysis
go_results_green <- enrichGO(
  gene         = entrez_ids_green,
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable     = TRUE
)

#No genes can be mapped

barplot(go_results_green, showCategory = 15, title = "Top GO Terms for Green Module")+
  theme(axis.text.y = element_text(size = 8))

dotplot(go_results_green, showCategory = 15, title = "Dotplot for Green Module GO Terms")+
  theme(axis.text.y = element_text(size = 8))

# Save gene set with Entrez IDs and gene symbols for Cytoscape
# Flatten the EntrezID list into a character string for the brown module
green_genes_with_entrez <- green_genes_with_entrez %>%
  mutate(EntrezID = sapply(EntrezID, function(x) paste(x, collapse = ",")))

# Save gene set with Entrez IDs and gene symbols for Cytoscape
write.csv(green_genes_with_entrez, file = "green_genes_with_entrez.csv", row.names = FALSE)

# Filter the genes for the "red" module
red_genes_with_entrez <- genes_by_module %>%
  filter(module1 == "red" & !is.na(EntrezID))

# Extract the Entrez IDs from the filtered genes (in this case, red module)
entrez_ids_red <- unlist(red_genes_with_entrez$EntrezID)

# Perform GO enrichment analysis
go_results_red <- enrichGO(
  gene         = entrez_ids_red,
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable     = TRUE
)

barplot(go_results_red, showCategory = 15, title = "Top GO Terms for Red Module")+
  theme(axis.text.y = element_text(size = 8))

dotplot(go_results_red, showCategory = 15, title = "Dotplot for Red Module GO Terms")+
  theme(axis.text.y = element_text(size = 8))

# Save gene set with Entrez IDs and gene symbols for Cytoscape
# Flatten the EntrezID list into a character string for the brown module
red_genes_with_entrez <- red_genes_with_entrez %>%
  mutate(EntrezID = sapply(EntrezID, function(x) paste(x, collapse = ",")))

# Save gene set with Entrez IDs and gene symbols for Cytoscape
write.csv(red_genes_with_entrez, file = "red_genes_with_entrez.csv", row.names = FALSE)

# Filter the genes for the "yellow" module
yellow_genes_with_entrez <- genes_by_module %>%
  filter(module1 == "yellow" & !is.na(EntrezID))

# Extract the Entrez IDs from the filtered genes (in this case, yellow module)
entrez_ids_yellow <- unlist(yellow_genes_with_entrez$EntrezID)

# Perform GO enrichment analysis
go_results_yellow <- enrichGO(
  gene         = entrez_ids_yellow,
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable     = TRUE
)

barplot(go_results_yellow, showCategory = 15, title = "Top GO Terms for Yellow Module")+
  theme(axis.text.y = element_text(size = 8))

dotplot(go_results_yellow, showCategory = 15, title = "Dotplot for Yellow Module GO Terms")+
  theme(axis.text.y = element_text(size = 8))

# Save gene set with Entrez IDs and gene symbols for Cytoscape
# Flatten the EntrezID list into a character string for the brown module
yellow_genes_with_entrez <- yellow_genes_with_entrez %>%
  mutate(EntrezID = sapply(EntrezID, function(x) paste(x, collapse = ",")))

# Save gene set with Entrez IDs and gene symbols for Cytoscape
write.csv(yellow_genes_with_entrez, file = "yellow_genes_with_entrez.csv", row.names = FALSE)


######
library(biomaRt)

# Use the appropriate mart for mouse
ensembl = useMart("ensembl", dataset = "mmusculus_gene_ensembl") # For mouse

# Load the data file with Ensembl IDs
edge_list <- read.csv("C:/Users/agsko/dev/pcm/WGCNA/edgelist_all.tsv", sep = "\t")

# Extract ENSEMBL IDs (assuming gene1 and gene2 contain the ENSEMBL IDs)
ensembl_ids <- unique(c(edge_list$gene1, edge_list$gene2))

# Perform the mapping from ENSEMBL to gene symbols
mapping <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
                 filters = 'ensembl_gene_id', 
                 values = ensembl_ids, 
                 mart = ensembl)

# Merge the mapping with the edgelist
edgelist_mapped <- merge(edge_list, mapping, by.x = 'gene1', by.y = 'ensembl_gene_id', all.x = TRUE)
edgelist_mapped <- merge(edgelist_mapped, mapping, by.x = 'gene2', by.y = 'ensembl_gene_id', all.x = TRUE, suffixes = c("_gene1", "_gene2"))

# Filter by turquoise and red modules if relevant
edgelist_filtered <- subset(edgelist_mapped, module1 == "turquoise" & module2 == "red")

# Save the new file
write.csv(edgelist_filtered, "C:/Users/agsko/dev/pcm/WGCNA/edgelist_with_gene_symbols.tsv", row.names = FALSE)

# Load your red genes CSV file
red_genes <- read.csv("C:/Users/agsko/dev/pcm/WGCNA/red_genes_with_entrez.csv")

# Check the structure of the 'red_genes' and 'all.gene.meta' objects
str(red_genes)
str(all.gene.meta)

merged_red_genes <- merge(red_genes, all.gene.meta, by.x = "gene1", by.y = "ID", all.x = TRUE)

# View the merged data to confirm the Gene.name column has been added
head(merged_red_genes)

# Save the updated dataset with the gene symbol added
write.csv(merged_red_genes, "C:/Users/agsko/dev/pcm/WGCNA/merged_red_genes_with_gene_symbol.csv", row.names = FALSE)

# Load your red genes CSV file
red_genes <- read.csv("C:/Users/agsko/dev/pcm/WGCNA/membership_turquoise_hub_genes.csv")