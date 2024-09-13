install.packages("goat")

library(goat)

# TODO: change the output directory 
# or set to getwd() to write output to the current working directory
output_dir = getwd()
setwd("C:/Users/agsko/dev/pcm/GOAT")
library(clusterProfiler)

# Specify the paths to your GMT files
m8_path <- "C:/Users/agsko/dev/pcm/m8.all.v2024.1.Mm.entrez.gmt"
mh_path <- "C:/Users/agsko/dev/pcm/mh.all.v2024.1.Mm.entrez.gmt"
m2_path <- "C:/Users/agsko/dev/pcm/m2.all.v2024.1.Mm.entrez.gmt"
m5_path <- "C:/Users/agsko/dev/pcm/m5.all.v2024.1.Mm.entrez.gmt"
m3_path <- "C:/Users/agsko/dev/pcm/m3.all.v2024.1.Mm.entrez.gmt"

# Load the GMT files using the appropriate function
m8_genesets <- read.gmt(m8_path)

m2_genesets <- read.gmt(m2_path)
m5_genesets <- read.gmt(m5_path)
m3_genesets <- read.gmt(m5_path)


library(dplyr)

# Function to transform a gene set data frame with correct data types
transform_geneset <- function(geneset_df, source_name) {
  geneset_df %>%
    group_by(term) %>%
    summarize(
      source = as.character(source_name),
      source_version = as.character("v2024.1"),  # Ensure this is a character
      id = as.character(unique(term)),
      name = as.character(unique(term)),
      genes = list(gene),
      ngenes = n()
    ) %>%
    ungroup()
}


# Transform all gene sets
mh_genesets_transformed <- transform_geneset(mh_genesets, "MH")
m5_genesets_transformed <- transform_geneset(m5_genesets, "M5")
m2_genesets_transformed <- transform_geneset(m2_genesets, "M2")
m8_genesets_transformed <- transform_geneset(m8_genesets, "M8")
m3_genesets_transformed <- transform_geneset(m3_genesets, "M3")

# Combine all transformed gene sets into a single data frame
combined_genesets_transformed <- bind_rows(
  mh_genesets_transformed,
  m5_genesets_transformed,
  m2_genesets_transformed,
  m3_genesets_transformed,
  m8_genesets_transformed
)

# Inspect the combined transformed data
str(combined_genesets_transformed)

getwd()

# List of file paths provided by the user
file_paths <- c(
  "C:/Users/agsko/dev/pcm/GOAT/Ure_DMXAAvsPBS_RT_DMXAA_mapped_cleaned_final.csv",
  "C:/Users/agsko/dev/pcm/GOAT/Ure_DMXAAvsUre_mapped_cleaned_final.csv",
  "C:/Users/agsko/dev/pcm/GOAT/Ure_DMXAAvsUre_RT_mapped_cleaned_final.csv",
  "C:/Users/agsko/dev/pcm/GOAT/Ure_RT_DMXAAvsPBS_RT_DMXAA_mapped_cleaned_final.csv",
  "C:/Users/agsko/dev/pcm/GOAT/Ure_RT_DMXAAvsUre_DMXAA_mapped_cleaned_final.csv",
  "C:/Users/agsko/dev/pcm/GOAT/Ure_RT_DMXAAvsUre_mapped_cleaned_final.csv",
  "C:/Users/agsko/dev/pcm/GOAT/Ure_RTvsPBS_RT_DMXAA_mapped_cleaned_final.csv",
  "C:/Users/agsko/dev/pcm/GOAT/Ure_RT_DMXAAvsUre_RT_mapped_cleaned_final.csv",
  "C:/Users/agsko/dev/pcm/GOAT/Ure_RTvsUre_mapped_cleaned_final.csv",
  "C:/Users/agsko/dev/pcm/GOAT/UrevsPBS_RT_DMXAA_mapped_cleaned_final.csv"
)

# Loop through each file, remove NA values, and save the cleaned file
for (file_path in file_paths) {
  # Read the CSV file
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Remove rows with NA values
  df_cleaned <- na.omit(df)
  
  # Save the cleaned data back to the same file
  write.csv(df_cleaned, file_path, row.names = FALSE)
}


#####################
# Define the list of file paths
file_paths <- c(
  "C:/Users/agsko/dev/pcm/UrevsPBS_RT_DMXAA.csv",
  "C:/Users/agsko/dev/pcm/Ure_DMXAAvsPBS_RT_DMXAA.csv",
  "C:/Users/agsko/dev/pcm/Ure_DMXAAvsUre.csv",
  "C:/Users/agsko/dev/pcm/Ure_DMXAAvsUre_RT.csv",
  "C:/Users/agsko/dev/pcm/Ure_RT_DMXAAvsPBS_RT_DMXAA.csv",
  "C:/Users/agsko/dev/pcm/Ure_RT_DMXAAvsUre.csv",
  "C:/Users/agsko/dev/pcm/Ure_RT_DMXAAvsUre_DMXAA.csv",
  "C:/Users/agsko/dev/pcm/Ure_RT_DMXAAvsUre_RT.csv",
  "C:/Users/agsko/dev/pcm/Ure_RTvsPBS_RT_DMXAA.csv",
  "C:/Users/agsko/dev/pcm/Ure_RTvsUre.csv"
)

# Function to clean data by removing duplicates and filtering out rows with effect size = 0
clean_data <- function(file_path) {
  # Load the dataset
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Remove duplicates, keeping the entry with the smallest p-value
  data <- data[order(data$gene, data$pvalue), ]  # Order by gene and pvalue
  data_unique <- data[!duplicated(data$gene), ]  # Remove duplicates, keeping the first (most significant)
  
  # Remove genes with effect size = 0
  data_unique <- data_unique[data_unique$effectsize != 0, ]
  
  # Check the number of genes remaining
  cat("Number of genes remaining after filtering in", file_path, ":", nrow(data_unique), "\n")
  
  # Save the cleaned data back to a new file
  cleaned_file_path <- gsub(".csv", "_cleaned.csv", file_path)
  write.csv(data_unique, cleaned_file_path, row.names = FALSE)
}

# Loop through each file and clean it
for (file_path in file_paths) {
  clean_data(file_path)
}

# Function to trim data down to 20,000 genes by removing those with the lowest absolute effect size and signif = FALSE
trim_data_to_20000 <- function(file_path) {
  # Load the cleaned dataset
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Calculate how many genes need to be removed
  num_genes_to_remove <- nrow(data) - 20000
  
  if (num_genes_to_remove > 0) {
    # Filter out non-significant genes
    non_signif_genes <- data[data$signif == "false", ]
    
    # Order by absolute effect size (lowest to highest)
    non_signif_genes <- non_signif_genes[order(abs(non_signif_genes$effectsize)), ]
    
    # Select genes to remove
    genes_to_remove <- non_signif_genes[1:num_genes_to_remove, ]
    
    # Remove these genes from the original data
    data <- data[!(data$gene %in% genes_to_remove$gene), ]
  }
  
  # Check the number of genes remaining
  cat("Number of genes remaining after trimming in", file_path, ":", nrow(data), "\n")
  
  # Save the trimmed data back to a new file
  trimmed_file_path <- gsub("_cleaned.csv", "_trimmed.csv", file_path)
  write.csv(data, trimmed_file_path, row.names = FALSE)
}

# List of all cleaned file paths
cleaned_file_paths <- c(
  "C:/Users/agsko/dev/pcm/UrevsPBS_RT_DMXAA_cleaned.csv",
  "C:/Users/agsko/dev/pcm/Ure_DMXAAvsPBS_RT_DMXAA_cleaned.csv",
  "C:/Users/agsko/dev/pcm/Ure_DMXAAvsUre_cleaned.csv",
  "C:/Users/agsko/dev/pcm/Ure_DMXAAvsUre_RT_cleaned.csv",
  "C:/Users/agsko/dev/pcm/Ure_RT_DMXAAvsPBS_RT_DMXAA_cleaned.csv",
  "C:/Users/agsko/dev/pcm/Ure_RT_DMXAAvsUre_cleaned.csv",
  "C:/Users/agsko/dev/pcm/Ure_RT_DMXAAvsUre_DMXAA_cleaned.csv",
  "C:/Users/agsko/dev/pcm/Ure_RT_DMXAAvsUre_RT_cleaned.csv",
  "C:/Users/agsko/dev/pcm/Ure_RTvsPBS_RT_DMXAA_cleaned.csv",
  "C:/Users/agsko/dev/pcm/Ure_RTvsUre_cleaned.csv"
)

# Loop through each file and trim it
for (file_path in cleaned_file_paths) {
  trim_data_to_20000(file_path)
}


#######################################

# Generate heatmaps for the clustered gene sets
plot_heatmap(clusters, output_dir)

# Generate lollipop plots for the reduced gene sets
plot_lollipop(result_reduced, output_dir, only_reduced = TRUE, plot_type = "lollipop", score_xaxis = "minlogp", score_color = "updown")

###############################################
# Load necessary libraries
library(goat)
library(dplyr)
library(openxlsx)

# Define the list of trimmed files
file_paths <- c(
  "C:/Users/agsko/dev/pcm/GOAT/Ure_DMXAAvsPBS_RT_DMXAA_mapped_cleaned_final.csv",
  "C:/Users/agsko/dev/pcm/GOAT/Ure_DMXAAvsUre_mapped_cleaned_final.csv",
  "C:/Users/agsko/dev/pcm/GOAT/Ure_DMXAAvsUre_RT_mapped_cleaned_final.csv",
  "C:/Users/agsko/dev/pcm/GOAT/Ure_RT_DMXAAvsPBS_RT_DMXAA_mapped_cleaned_final.csv",
  "C:/Users/agsko/dev/pcm/GOAT/Ure_RT_DMXAAvsUre_DMXAA_mapped_cleaned_final.csv",
  "C:/Users/agsko/dev/pcm/GOAT/Ure_RT_DMXAAvsUre_mapped_cleaned_final.csv",
  "C:/Users/agsko/dev/pcm/GOAT/Ure_RTvsPBS_RT_DMXAA_mapped_cleaned_final.csv",
  "C:/Users/agsko/dev/pcm/GOAT/Ure_RT_DMXAAvsUre_RT_mapped_cleaned_final.csv",
  "C:/Users/agsko/dev/pcm/GOAT/Ure_RTvsUre_mapped_cleaned_final.csv",
  "C:/Users/agsko/dev/pcm/GOAT/UrevsPBS_RT_DMXAA_mapped_cleaned_final.csv"
)

# Function to perform the GOAT analysis
perform_goat_analysis <- function(file_path) {
  # Load the dataset
  datasets <- read.csv(file_path, row.names = 1)
  
  # Create the genelist by extracting relevant columns
  genelist <- data.frame(
    gene = rownames(datasets),          # Extract the row names as 'gene'
    effectsize = datasets$effectsize,   # Extract the 'effectsize' column
    signif = as.logical(datasets$signif)  # Convert 'signif' to logical (TRUE/FALSE)
  )
  
  # Create a folder for the contrast
  contrast_name <- sub("_trimmed.csv", "", basename(file_path))
  contrast_folder <- file.path("C:/Users/agsko/dev/pcm/GOAT", contrast_name)
  dir.create(contrast_folder, showWarnings = FALSE)
  
  # Filter gene sets
  genesets_filtered <- filter_genesets(combined_genesets_transformed, genelist, min_overlap = 10, max_overlap = 1500)
  
  # Run the GOAT analysis
  result <- test_genesets(genesets_filtered, genelist, method = "goat", score_type = "effectsize", padj_method = "fdr", padj_cutoff = 0.05)
  
  # Print the significant GO term counts, per geneset 'source' (CC/BP/MF), to console
  print(result |> group_by(source) |> summarise(signif_count = sum(signif), .groups="drop"))
  
  # Filter the results for significant gene sets and select relevant columns
  significant_results <- result %>% 
    filter(signif) %>% 
    dplyr::select(source, name, ngenes, pvalue_adjust)
  
  # Save the significant results to a CSV file
  output_filename <- file.path(contrast_folder, paste0(contrast_name, "_GOAT.csv"))
  write.csv(significant_results, output_filename, row.names = FALSE)
  
  # Get the top 25 significant results
  top_25_results <- head(significant_results, 25)
  
  # Print and save the top 25 significant results
  top_25_filename <- file.path(contrast_folder, paste0(contrast_name, "_GOAT_25.csv"))
  write.csv(top_25_results, top_25_filename, row.names = FALSE)
  
  # Print significant gene set counts before and after simplification
  print(result |> dplyr::filter(signif) |> dplyr::count(source))
  print(result_reduced |> dplyr::filter(signif_and_reduced) |> dplyr::count(source))
  

# Loop through each trimmed file and perform the GOAT analysis
for (file_path in trimmed_files) {
  perform_goat_analysis(file_path)
}
  
#########################

#Ure_DMXAAvsUre
  
  # Run the GOAT analysis
  result_Ure_DMXAAvsUre <- test_genesets(genesets_filtered, genelist, method = "goat", score_type = "effectsize", padj_method = "fdr", padj_cutoff = 0.05)
  
  # Generate clusters
  cluster_Ure_DMXAAvsUre <- cluster_genesets(result_Ure_DMXAAvsUre, genelist)
  
  # Reduce gene sets to remove overlap
  result_reduced_Ure_DMXAAvsUre <- reduce_genesets(cluster_Ure_DMXAAvsUre, simscore_threshold = 0.9, universe_fraction = 0.25, signifgenes_fraction = 0.9)
  
  # Specify the output directory
  output_dir <- "C:/Users/agsko/dev/pcm/Ure_DMXAAvsUre"
  
  # Generate and save the heatmap plot
  plot_heatmap(cluster_Ure_DMXAAvsUre, output_dir)
  
  # Generate and save the lollipop plot
  plot_lollipop(result_reduced_Ure_DMXAAvsUre, output_dir, only_reduced = TRUE, plot_type = "lollipop", score_xaxis = "minlogp", score_color = "updown")
  
  # Ure_RT vs Ure
  
  # Run the GOAT analysis
  result_Ure_RTvsUre <- test_genesets(genesets_filtered, genelist, method = "goat", score_type = "effectsize", padj_method = "fdr", padj_cutoff = 0.05)
  
  # Generate clusters
  cluster_Ure_RTvsUre <- cluster_genesets(result_Ure_RTvsUre, genelist)
  
  # Reduce gene sets to remove overlap
  result_reduced_Ure_RTvsUre <- reduce_genesets(cluster_Ure_RTvsUre, simscore_threshold = 0.9, universe_fraction = 0.25, signifgenes_fraction = 0.9)
  
  # Specify the output directory
  output_dir <- "C:/Users/agsko/dev/pcm/Ure_RTvsUre"
  
  # Generate and save the heatmap plot
  plot_heatmap(cluster_Ure_RTvsUre, output_dir)
  
  # Generate and save the lollipop plot
  plot_lollipop(result_reduced_Ure_RTvsUre, output_dir, only_reduced = TRUE, plot_type = "lollipop", score_xaxis = "minlogp", score_color = "updown")
  
  # Ure_RT_DMXAAvsUre_RT
  # Run the GOAT analysis
  result_Ure_RT_DMXAAvsUre_RT <- test_genesets(genesets_filtered, genelist, method = "goat", score_type = "effectsize", padj_method = "fdr", padj_cutoff = 0.05)
  
  # Generate clusters
  cluster_Ure_RT_DMXAAvsUre_RT <- cluster_genesets(result_Ure_RT_DMXAAvsUre_RT, genelist)
  
  # Reduce gene sets to remove overlap
  result_reduced_Ure_RT_DMXAAvsUre_RT <- reduce_genesets(cluster_Ure_RT_DMXAAvsUre_RT, simscore_threshold = 0.9, universe_fraction = 0.25, signifgenes_fraction = 0.9)
  
  # Specify the output directory
  output_dir <- "C:/Users/agsko/dev/pcm/Ure_RT_DMXAAvsUre_RT"
  
  # Generate and save the heatmap plot
  plot_heatmap(cluster_Ure_RT_DMXAAvsUre_RT, output_dir)
  
  # Generate and save the lollipop plot
  plot_lollipop(result_reduced_Ure_RT_DMXAAvsUre_RT, output_dir, only_reduced = TRUE, plot_type = "lollipop", score_xaxis = "minlogp", score_color = "updown")
  
  #Ure_RT_DMXAAvsUre
  # Run the GOAT analysis
  result_Ure_RT_DMXAAvsUre <- test_genesets(genesets_filtered, genelist, method = "goat", score_type = "effectsize", padj_method = "fdr", padj_cutoff = 0.05)
  
  # Generate clusters
  cluster_Ure_RT_DMXAAvsUre <- cluster_genesets(result_Ure_RT_DMXAAvsUre, genelist)
  
  # Reduce gene sets to remove overlap
  result_reduced_Ure_RT_DMXAAvsUre <- reduce_genesets(cluster_Ure_RT_DMXAAvsUre, simscore_threshold = 0.9, universe_fraction = 0.25, signifgenes_fraction = 0.9)
  
  # Specify the output directory
  output_dir <- "C:/Users/agsko/dev/pcm/Ure_RT_DMXAAvsUre"
  
  # Generate and save the heatmap plot
  plot_heatmap(cluster_Ure_RT_DMXAAvsUre, output_dir)
  
  # Generate and save the lollipop plot
  plot_lollipop(result_reduced_Ure_RT_DMXAAvsUre, output_dir, only_reduced = TRUE, plot_type = "lollipop", score_xaxis = "minlogp", score_color = "updown")
  
  #UrevsPBS_RT_DMXAA
  # Run the GOAT analysis
  result_UrevsPBS_RT_DMXAA <- test_genesets(genesets_filtered, genelist, method = "goat", score_type = "effectsize", padj_method = "fdr", padj_cutoff = 0.05)
  
  # Generate clusters
  cluster_UrevsPBS_RT_DMXAA <- cluster_genesets(result_UrevsPBS_RT_DMXAA, genelist)
  
  # Reduce gene sets to remove overlap
  result_reduced_UrevsPBS_RT_DMXAA <- reduce_genesets(cluster_UrevsPBS_RT_DMXAA, simscore_threshold = 0.9, universe_fraction = 0.25, signifgenes_fraction = 0.9)
  
  # Specify the output directory
  output_dir <- "C:/Users/agsko/dev/pcm/UrevsPBS_RT_DMXAA"
  
  # Generate and save the heatmap plot
  plot_heatmap(cluster_UrevsPBS_RT_DMXAA, output_dir)
  
  # Generate and save the lollipop plot
  plot_lollipop(result_reduced_UrevsPBS_RT_DMXAA, output_dir, only_reduced = TRUE, plot_type = "lollipop", score_xaxis = "minlogp", score_color = "updown")
  
  # Ure_DMXAAvsPBS_RT_DMXAA
  # Run the GOAT analysis
  result_Ure_DMXAAvsPBS_RT_DMXAA <- test_genesets(genesets_filtered, genelist, method = "goat", score_type = "effectsize", padj_method = "fdr", padj_cutoff = 0.05)
  
  # Generate clusters
  cluster_Ure_DMXAAvsPBS_RT_DMXAA <- cluster_genesets(result_Ure_DMXAAvsPBS_RT_DMXAA, genelist)
  
  # Reduce gene sets to remove overlap
  result_reduced_Ure_DMXAAvsPBS_RT_DMXAA <- reduce_genesets(cluster_Ure_DMXAAvsPBS_RT_DMXAA, simscore_threshold = 0.9, universe_fraction = 0.25, signifgenes_fraction = 0.9)
  
  # Specify the output directory
  output_dir <- "C:/Users/agsko/dev/pcm/Ure_DMXAAvsPBS_RT_DMXAA"
  
  # Generate and save the heatmap plot
  plot_heatmap(cluster_Ure_DMXAAvsPBS_RT_DMXAA, output_dir)
  
  # Generate and save the lollipop plot
  plot_lollipop(result_reduced_Ure_DMXAAvsPBS_RT_DMXAA, output_dir, only_reduced = TRUE, plot_type = "lollipop", score_xaxis = "minlogp", score_color = "updown")
  
  #Ure_RTvsPBS_RT_DMXAA
  # Run the GOAT analysis
  result_Ure_RTvsPBS_RT_DMXAA <- test_genesets(genesets_filtered, genelist, method = "goat", score_type = "effectsize", padj_method = "fdr", padj_cutoff = 0.05)
  
  # Generate clusters
  cluster_Ure_RTvsPBS_RT_DMXAA <- cluster_genesets(result_Ure_RTvsPBS_RT_DMXAA, genelist)
  
  # Reduce gene sets to remove overlap
  result_reduced_Ure_RTvsPBS_RT_DMXAA <- reduce_genesets(cluster_Ure_RTvsPBS_RT_DMXAA, simscore_threshold = 0.9, universe_fraction = 0.25, signifgenes_fraction = 0.9)
  
  # Specify the output directory
  output_dir <- "C:/Users/agsko/dev/pcm/Ure_RTvsPBS_RT_DMXAA"
  
  # Generate and save the heatmap plot
  plot_heatmap(cluster_Ure_RTvsPBS_RT_DMXAA, output_dir)
  
  # Generate and save the lollipop plot
  plot_lollipop(result_reduced_Ure_RTvsPBS_RT_DMXAA, output_dir, only_reduced = TRUE, plot_type = "lollipop", score_xaxis = "minlogp", score_color = "updown")
  
  #Ure_RT_DMXAAvsUre_RT
  # Run the GOAT analysis
  result_Ure_RT_DMXAAvsUre_RT <- test_genesets(genesets_filtered, genelist, method = "goat", score_type = "effectsize", padj_method = "fdr", padj_cutoff = 0.05)
  
  # Generate clusters
  cluster_Ure_RT_DMXAAvsUre_RT <- cluster_genesets(result_Ure_RT_DMXAAvsUre_RT, genelist)
  
  # Reduce gene sets to remove overlap
  result_reduced_Ure_RT_DMXAAvsUre_RT <- reduce_genesets(cluster_Ure_RT_DMXAAvsUre_RT, simscore_threshold = 0.9, universe_fraction = 0.25, signifgenes_fraction = 0.9)
  
  # Specify the output directory
  output_dir <- "C:/Users/agsko/dev/pcm/Ure_RT_DMXAAvsUre_RT"
  
  # Generate and save the heatmap plot
  plot_heatmap(cluster_Ure_RT_DMXAAvsUre_RT, output_dir)
  
  # Generate and save the lollipop plot
  plot_lollipop(result_reduced_Ure_RT_DMXAAvsUre_RT, output_dir, only_reduced = TRUE, plot_type = "lollipop", score_xaxis = "minlogp",
                
  #Ure_RT_DMXAAvsPBS_RT_DMXAA
  # Run the GOAT analysis
  result_Ure_RT_DMXAAvsPBS_RT_DMXAA <- test_genesets(genesets_filtered, genelist, method = "goat", score_type = "effectsize", padj_method = "fdr", padj_cutoff = 0.05)
  
  # Generate clusters
  cluster_Ure_RT_DMXAAvsPBS_RT_DMXAA <- cluster_genesets(result_Ure_RT_DMXAAvsPBS_RT_DMXAA, genelist)
  
  # Reduce gene sets to remove overlap
  result_reduced_Ure_RT_DMXAAvsPBS_RT_DMXAA <- reduce_genesets(cluster_Ure_RT_DMXAAvsPBS_RT_DMXAA, simscore_threshold = 0.9, universe_fraction = 0.25, signifgenes_fraction = 0.9)
  
  # Specify the output directory
  output_dir <- "C:/Users/agsko/dev/pcm/Ure_RT_DMXAAvsPBS_RT_DMXAA"
  
  # Generate and save the heatmap plot
  plot_heatmap(cluster_Ure_RT_DMXAAvsPBS_RT_DMXAA, output_dir)
  
  # Generate and save the lollipop plot
  plot_lollipop(result_reduced_Ure_RT_DMXAAvsPBS_RT_DMXAA, output_dir, only_reduced = TRUE, plot_type = "lollipop", score_xaxis = "minlogp", score_color = "updown")
  
  #Ure_DMXAAvsUre_RT
  # Run the GOAT analysis
  result_Ure_DMXAAvsUre_RT <- test_genesets(genesets_filtered, genelist, method = "goat", score_type = "effectsize", padj_method = "fdr", padj_cutoff = 0.05)
  
  # Generate clusters
  cluster_Ure_DMXAAvsUre_RT <- cluster_genesets(result_Ure_DMXAAvsUre_RT, genelist)
  
  # Reduce gene sets to remove overlap
  result_reduced_Ure_DMXAAvsUre_RT <- reduce_genesets(cluster_Ure_DMXAAvsUre_RT, simscore_threshold = 0.9, universe_fraction = 0.25, signifgenes_fraction = 0.9)
  
  # Specify the output directory
  output_dir <- "C:/Users/agsko/dev/pcm/Ure_DMXAAvsUre_RT"
  
  # Generate and save the heatmap plot
  plot_heatmap(cluster_Ure_DMXAAvsUre_RT, output_dir)
  
  # Generate and save the lollipop plot
  plot_lollipop(result_reduced_Ure_DMXAAvsUre_RT, output_dir, only_reduced = TRUE, plot_type = "lollipop", score_xaxis = "minlogp", score_color = "updown")
  
  #Ure_RT_DMXAAvsUre_DMXAA
  # Run the GOAT analysis
  result_Ure_RT_DMXAAvsUre_DMXAA <- test_genesets(genesets_filtered, genelist, method = "goat", score_type = "effectsize", padj_method = "fdr", padj_cutoff = 0.05)
  
  # Generate clusters
  cluster_Ure_RT_DMXAAvsUre_DMXAA <- cluster_genesets(result_Ure_RT_DMXAAvsUre_DMXAA, genelist)
  
  # Reduce gene sets to remove overlap
  result_reduced_Ure_RT_DMXAAvsUre_DMXAA <- reduce_genesets(cluster_Ure_RT_DMXAAvsUre_DMXAA, simscore_threshold = 0.9, universe_fraction = 0.25, signifgenes_fraction = 0.9)
  
  # Specify the output directory
  output_dir <- "C:/Users/agsko/dev/pcm/Ure_RT_DMXAAvsUre_DMXAA"
  
  # Generate and save the heatmap plot
  plot_heatmap(cluster_Ure_RT_DMXAAvsUre_DMXAA, output_dir)
  
  # Generate and save the lollipop plot
  plot_lollipop(result_reduced_Ure_RT_DMXAAvsUre_DMXAA, output_dir, only_reduced = TRUE, plot_type = "lollipop", score_xaxis = "minlogp", score_color = "updown")
  