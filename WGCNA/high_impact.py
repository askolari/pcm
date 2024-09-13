import pandas as pd

# File paths
high_impact_file = r"/mnt/c/Users/agsko/dev/pcm/WES_explor/high_impact_var_ure.xlsx"
upregulated_genes_file = r"/mnt/c/Users/agsko/dev/pcm/predictions/merged_upregulated_genes.csv"
downregulated_genes_file = r"/mnt/c/Users/agsko/dev/pcm/predictions/merged_downregulated_genes.csv"

# Output file paths
output_upregulated = r"/mnt/c/Users/agsko/dev/pcm/WGCNA/high_impact_upregulated_genes.csv"
output_downregulated = r"/mnt/c/Users/agsko/dev/pcm/WGCNA/high_impact_downregulated_genes.csv"

# Read the high-impact variants Excel file
high_impact_df = pd.read_excel(high_impact_file)

# Assuming 'Gene.name' is the column containing the gene names
high_impact_genes = high_impact_df['Gene.name'].unique()

# Read the upregulated and downregulated gene CSV files
upregulated_genes_df = pd.read_csv(upregulated_genes_file)
downregulated_genes_df = pd.read_csv(downregulated_genes_file)

# Filter for common occurrences in upregulated genes
upregulated_common = upregulated_genes_df[upregulated_genes_df['Gene.name'].isin(high_impact_genes)]

# Filter for common occurrences in downregulated genes
downregulated_common = downregulated_genes_df[downregulated_genes_df['Gene.name'].isin(high_impact_genes)]

# Save the filtered data to new CSV files
upregulated_common.to_csv(output_upregulated, index=False)
downregulated_common.to_csv(output_downregulated, index=False)

print("Filtered upregulated and downregulated genes saved.")
