# List of file paths
file_paths <- c(
  "C:/Users/agsko/dev/pcm/EdgeR/Venn and common genes/common_ure_ure_rt_genes_upregulated_with_names.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/Venn and common genes/common_all_ure_ure_rt_ure_dmxaa_upregulated_with_names.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/Venn and common genes/common_ure_dmxaa_ure_rt_dmxaa_genes_upregulated_with_names.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/Venn and common genes/common_ure_rt_ure_dmxaa_genes_upregulated_with_names.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/Venn and common genes/common_ure_rt_ure_rt_dmxaa_genes_upregulated_with_names.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/Venn and common genes/common_ure_ure_dmxaa_genes_upregulated_with_names.csv"
)

# Create a vector for the group names based on file names
group_names <- c(
  "ure_ure_rt",
  "all_ure_ure_rt_ure_dmxaa",
  "ure_dmxaa_ure_rt_dmxaa",
  "ure_rt_ure_dmxaa",
  "ure_rt_ure_rt_dmxaa",
  "ure_ure_dmxaa"
)

# Load and merge the CSV files, adding a new column for group
merged_data <- lapply(seq_along(file_paths), function(i) {
  data <- read.csv(file_paths[i])
  data$Group <- group_names[i]
  return(data)
})

# Combine all data frames into one
merged_data_df <- do.call(rbind, merged_data)

# Save the combined data into a new CSV file
write.csv(merged_data_df, file = "merged_common_genes_upregulated_with_groups.csv", row.names = FALSE)

# Preview the combined data
head(merged_data_df)
