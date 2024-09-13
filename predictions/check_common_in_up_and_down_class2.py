import pandas as pd

# Define the file paths
common_peptides_file = "/mnt/c/Users/agsko/dev/pcm/predictions/common_peptides_in_different_samples_class2.csv"
upregulated_genes_file = "/mnt/c/Users/agsko/dev/pcm/predictions/merged_upregulated_genes.csv"
downregulated_genes_file = "/mnt/c/Users/agsko/dev/pcm/predictions/merged_downregulated_genes.csv"

# Load the data
common_peptides_df = pd.read_csv(common_peptides_file)
upregulated_genes_df = pd.read_csv(upregulated_genes_file)
downregulated_genes_df = pd.read_csv(downregulated_genes_file)

# Check if any genes from the common peptides file are present in either the upregulated or downregulated genes
common_genes_upregulated = pd.merge(common_peptides_df, upregulated_genes_df, on="Gene.name", how="inner")
common_genes_downregulated = pd.merge(common_peptides_df, downregulated_genes_df, on="Gene.name", how="inner")

# Concatenate the results from both upregulated and downregulated matches
common_genes = pd.concat([common_genes_upregulated, common_genes_downregulated], ignore_index=True)

# Save the result to a new file
output_file = "/mnt/c/Users/agsko/dev/pcm/predictions/common_genes_in_peptides_upregulated_or_downregulated_class2.csv"
common_genes.to_csv(output_file, index=False)

print(f"Common genes saved to {output_file}")

# Check if any genes from the common peptides file are present in the upregulated genes
common_genes_upregulated = pd.merge(common_peptides_df, upregulated_genes_df, on="Gene.name", how="inner")

# Save the result to a new file
output_file_upregulated = "/mnt/c/Users/agsko/dev/pcm/predictions/common_genes_in_peptides_upregulated_class2.csv"
common_genes_upregulated.to_csv(output_file_upregulated, index=False)

print(f"Common genes saved to {output_file_upregulated}")