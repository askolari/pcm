import pandas as pd

# File paths
cytoscape_file = r"/mnt/c/Users/agsko/dev/pcm/WGCNA/cytoscape_node_list.csv"
turquoise_genes_file = r"/mnt/c/Users/agsko/dev/pcm/WGCNA/turquoise_genes_with_entrez.txt"

# Output file path
output_file = r"/mnt/c/Users/agsko/dev/pcm/WGCNA/turquoise_genes_with_cytoscape_data.csv"

# Read the files into pandas dataframes
cytoscape_df = pd.read_csv(cytoscape_file)
turquoise_genes_df = pd.read_csv(turquoise_genes_file, sep='\t')

# Find the common occurrences in 'gene1' and 'gene' columns
merged_df = pd.merge(turquoise_genes_df, cytoscape_df, left_on='gene1', right_on='gene1', how='inner')

# Save the new dataframe with merged columns to a new file
merged_df.to_csv(output_file, index=False)

print(f"Matched and merged data saved to: {output_file}")
