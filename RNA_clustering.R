# Perform voom transformation
v <- voom(dge_filtered, design, plot = TRUE)

# Extract the voom-transformed expression values (logCPM)
voom_logcpm <- v$E

# Step 1: Perform k-means clustering without centering the logCPM data
scree_data_no_centering <- data.frame(k = 1:20)

km_list_no_centering <- bpmapply(kmeans, centers = scree_data_no_centering$k,
                                 MoreArgs = list(x = voom_logcpm, nstart = 5),
                                 SIMPLIFY = FALSE)
scree_data_no_centering$wss <- sapply(km_list_no_centering, function(km) km$tot.withinss)

# Step 2: Create the scree plot without centering
ggplot(scree_data_no_centering) + 
  aes(x = k, y = wss) +
  geom_line(alpha = 0.5) + 
  geom_point() +
  ggtitle("K-means Scree Plot (No Centering)") +
  labs(x = "Number of Clusters (k)", y = "Total Within-Cluster Sum of Squares") +
  theme_minimal()

# Step 1: Center the logCPM data by subtracting the row means
logCPM_centered <- t(scale(t(voom_logcpm), scale = FALSE))

# Step 2: Perform k-means clustering with centering the logCPM data
scree_data_centering <- data.frame(k = 1:20)

km_list_centering <- bpmapply(kmeans, centers = scree_data_centering$k,
                              MoreArgs = list(x = logCPM_centered, nstart = 5),
                              SIMPLIFY = FALSE)
scree_data_centering$wss <- sapply(km_list_centering, function(km) km$tot.withinss)

# Step 3: Create the scree plot with centering
ggplot(scree_data_centering) + 
  aes(x = k, y = wss) +
  geom_line(alpha = 0.5) + 
  geom_point() +
  ggtitle("K-means Scree Plot (With Centering)") +
  labs(x = "Number of Clusters (k)", y = "Total Within-Cluster Sum of Squares") +
  theme_minimal()

#We will continue with hierarchical clustering and 3 clusters, and then 4

# Set the number of clusters to 3
set.seed(123)  # Set a seed for reproducibility
kmeans_result_3 <- kmeans(voom_logcpm, centers = 3, nstart = 5)

# View cluster assignments
table(kmeans_result_3$cluster)

# Plot the clusters (PCA for visualization)
pca_result <- prcomp(voom_logcpm, scale. = TRUE)
pca_data <- data.frame(PC1 = pca_result$x[,1], PC2 = pca_result$x[,2], 
                       Cluster = as.factor(kmeans_result_3$cluster))
ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 3) +
  labs(title = "K-means Clustering with 3 Clusters (PCA plot)") +
  theme_minimal()

# Set the number of clusters to 4
set.seed(123)  # Set a seed for reproducibility
kmeans_result_4 <- kmeans(voom_logcpm, centers = 4, nstart = 5)

# View cluster assignments
table(kmeans_result_4$cluster)

# Plot the clusters (PCA for visualization)
pca_data <- data.frame(PC1 = pca_result$x[,1], PC2 = pca_result$x[,2], 
                       Cluster = as.factor(kmeans_result_4$cluster))
ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 3) +
  labs(title = "K-means Clustering with 4 Clusters (PCA plot)") +
  theme_minimal()

# Transpose voom_logcpm so that samples are clustered instead of genes
voom_logcpm_t <- t(voom_logcpm)

# Perform k-means clustering on samples (transposed matrix)
kmeans_result_samples_4 <- kmeans(voom_logcpm_t, centers = 4, nstart = 5)

# Now you can create the contingency table with clusters and conditions
contingency_table <- table(kmeans_result_samples_4$cluster, metadata_combined$Condition)
print(contingency_table)

# Perform PCA on the transposed voom_logcpm matrix (clustering was done on samples)
pca_result <- prcomp(voom_logcpm_t, scale. = TRUE)

# Create a data frame for plotting PCA results
pca_data <- as.data.frame(pca_result$x)
pca_data$Condition <- metadata_combined$Condition

# Plot PCA with samples colored by Condition
ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot Colored by Condition") +
  theme_minimal()

# Perform a Chi-squared test to check for association between clusters and condition
chi_sq_test <- chisq.test(contingency_table)
print(chi_sq_test)

#Pearson's Chi-squared test
#data:  contingency_table
#X-squared = 21, df = 12, p-value = 0.05038

# Load necessary libraries
library(cluster)
install.packages("factoextra")
library(factoextra)

# 1. Within-cluster sum of squares (WSS) for the 4-cluster model
wss_4_clusters <- sum(kmeans_result_4$withinss)
#WSS: 196118.3 

# 2. Silhouette score for the 4-cluster model
silhouette_4_clusters <- silhouette(kmeans_result_4$cluster, dist(voom_logcpm))
avg_silhouette_4 <- mean(silhouette_4_clusters[, 3])
#Average Silhouette Score: 0.4857938 

# Install clusterSim if not already installed
if (!require(clusterSim)) {
  install.packages("clusterSim")
}

# Load the clusterSim package
library(clusterSim)

# Calculate Davies-Bouldin Index for the 4-cluster model
db_index_4 <- index.DB(voom_logcpm, kmeans_result_4$cluster)$DB
cat("Davies-Bouldin Index for 4 clusters:", db_index_4, "\n")
#Davies-Bouldin Index for 4 clusters: 0.7197708 


install.packages("mclust")
library(mclust)

# Assuming you have the true labels in metadata_combined$Condition
# Load the necessary library
library(cluster)

# Compute the Adjusted Rand Index for sample clustering
ari_samples_4 <- adjustedRandIndex(kmeans_result_samples_4$cluster, metadata_combined$Condition)

# Print the ARI value
cat("Adjusted Rand Index for sample clustering (4 clusters):", ari_samples_4, "\n")
#Adjusted Rand Index for sample clustering (4 clusters): 0.2285714 

# Print metrics for the 4-cluster model
cat("WSS:", wss_4_clusters, "\n")
cat("Average Silhouette Score:", avg_silhouette_4, "\n")

#Hierarchical clustering

# Calculate the variance for each gene
rv <- apply(voom_logcpm, 1, var)

# Select the top 5000 most variable genes
topvar <- head(order(rv, decreasing = TRUE), 5000)

# Center the logCPM data for the top 5000 variable genes
logCPM_centered <- voom_logcpm[topvar, ] - rowMeans(voom_logcpm[topvar, ])

library(pheatmap)

# Create a color palette
heat_pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100)

# Perform hierarchical clustering and visualize it as a heatmap
pheatmap(logCPM_centered,
         annotation_col = metadata_combined[c("Condition", "RT_status")],
         show_rownames = FALSE,
         show_colnames = FALSE,
         cluster_cols = TRUE,  # Cluster samples (columns)
         clustering_distance_rows = "correlation",  # Cluster genes (rows) based on correlation
         color = heat_pal,
         breaks = seq(-6, 6, length.out = length(heat_pal) + 1))

# Calculate the variance for each gene
gene_variances <- apply(voom_logcpm, 1, var)

# Select the top 500 most variable genes
top_500_genes <- names(sort(gene_variances, decreasing = TRUE))[1:500]

# Subset the voom_logcpm matrix to include only the top 500 most variable genes
voom_logcpm_top500 <- voom_logcpm[top_500_genes, ]

# Subtract the row means using basic matrix operations
voom_logcpm_centered <- voom_logcpm_top500 - rowMeans(voom_logcpm_top500)

# Plot the simplified hierarchical clustering heatmap with Sample ID annotations
pheatmap(voom_logcpm_centered,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         show_rownames = FALSE, show_colnames = FALSE,
         cluster_cols = FALSE,
         clustering_distance_rows = "correlation",
         color = heat_pal, breaks = seq(-6, 6, length.out = length(heat_pal) + 1))

# Prepare annotation data for the columns (samples)
annotation_col <- metadata_combined[, c("Condition", "RT_status")]

# Define the colors for the annotations
library(RColorBrewer)
annotation_colors <- list(
  Condition = brewer.pal(length(unique(metadata_combined$Condition)), "Set1"),
  RT_status = c("RT" = "#4daf4a", "No_RT" = "#984ea3")  # Custom colors for RT status
)

# Assign names to annotation colors to match levels of the factors
names(annotation_colors$Condition) <- levels(metadata_combined$Condition)

# Plot the simplified hierarchical clustering heatmap with Condition and RT_status annotations
pheatmap(voom_logcpm_centered,
         annotation_col = annotation_col,  # Add the Condition and RT_status annotations
         annotation_colors = annotation_colors,  # Use the defined colors for annotations
         show_rownames = FALSE, show_colnames = FALSE,  # Hide row and column names
         cluster_cols = TRUE,  # Cluster the columns (samples)
         clustering_distance_rows = "correlation",  # Cluster the rows (genes) by correlation
         color = heat_pal,  # Use the color palette defined earlier
         breaks = seq(-6, 6, length.out = length(heat_pal) + 1))  # Set breaks for the color scale


# Load necessary libraries
library(ggplot2)
library(dplyr)

# Load the data
file_path <- "C:/Users/agsko/dev/pcm/predictions/combined_class1_predictions.csv"
data <- read.csv(file_path)

# Assign conditions based on Sample_ID
data$Condition <- case_when(
  data$Sample_ID == "S1" ~ "control",
  data$Sample_ID == "S2" ~ "control",
  data$Sample_ID == "S3" ~ "ure_mock_IR",
  data$Sample_ID == "S4" ~ "ure_mock_IR",
  data$Sample_ID == "S5" ~ "ure_mock_IR",
  data$Sample_ID == "S6" ~ "ure_mock_IR",
  data$Sample_ID == "S7" ~ "ure_DMXAA",
  data$Sample_ID == "S8" ~ "ure_DMXAA",
  data$Sample_ID == "S9" ~ "ure_DMXAA",
  data$Sample_ID == "S10" ~ "ure_DMXAA",
  data$Sample_ID == "S11" ~ "Ure_IR",
  data$Sample_ID == "S12" ~ "Ure_IR",
  data$Sample_ID == "S13" ~ "Ure_IR",
  data$Sample_ID == "S14" ~ "Ure_IR",
  data$Sample_ID == "S15" ~ "PBS_IR_DMXAA"
)

# Count the number of peptides predicted by Sample_ID and Condition
peptide_counts <- data %>%
  group_by(Sample_ID, Condition) %>%
  summarise(n_peptides = n())

# Create the bar plot
ggplot(peptide_counts, aes(x = Sample_ID, y = n_peptides, fill = Condition)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Number of Class 1 Predicted Peptides\nby Sample and Condition",
       x = "Sample ID", 
       y = "Number of Predicted Peptides") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2") # You can customize the color palette if needed

# Save the dataset with the new 'Condition' column
output_file_with_conditions <- "C:/Users/agsko/dev/pcm/predictions/combined_class1_predictions_with_conditions.csv"
write.csv(data, output_file_with_conditions, row.names = FALSE)

# Confirm the file was saved
print(paste("Dataset with 'Condition' column saved to:", output_file_with_conditions))

#Repeat for class2, class1_IC500 and class2_IC500