import pandas as pd

# Define file paths
matched_features_class1_file = "/mnt/c/Users/agsko/dev/pcm/predictions/matched_features_result.csv"
matched_features_class2_file = "/mnt/c/Users/agsko/dev/pcm/predictions/matched_features_class2_result.csv"
merged_upregulated_genes_file = "/mnt/c/Users/agsko/dev/pcm/predictions/merged_upregulated_genes.csv"

# Load the data
matched_features_class1_df = pd.read_csv(matched_features_class1_file)
matched_features_class2_df = pd.read_csv(matched_features_class2_file)
merged_upregulated_genes_df = pd.read_csv(merged_upregulated_genes_file)

# Merge matched features with upregulated genes based on the 'Gene' column in the first files and 'Gene.name' in the merged_upregulated_genes.csv file
prioritised_peptides_class1 = pd.merge(matched_features_class1_df, merged_upregulated_genes_df, left_on='Gene', right_on='Gene.name', how='inner')
prioritised_peptides_class2 = pd.merge(matched_features_class2_df, merged_upregulated_genes_df, left_on='Gene', right_on='Gene.name', how='inner')

# Save the results to new files
output_class1_file = "/mnt/c/Users/agsko/dev/pcm/predictions/prioritised_peptides_by_upreg_class1.csv"
output_class2_file = "/mnt/c/Users/agsko/dev/pcm/predictions/prioritised_peptides_by_upreg_class2.csv"

prioritised_peptides_class1.to_csv(output_class1_file, index=False)
prioritised_peptides_class2.to_csv(output_class2_file, index=False)

print(f"Prioritised peptides for class 1 saved to {output_class1_file}")
print(f"Prioritised peptides for class 2 saved to {output_class2_file}")