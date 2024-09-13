import pandas as pd

# Define file paths
edgelist_path = '/mnt/c/Users/agsko/dev/pcm/WGCNA/edgelist_all.tsv'
mapping_path = '/mnt/c/Users/agsko/dev/pcm/WGCNA/ensembl_to_gene_symbol.tsv'

# Load edgelist file
edgelist_df = pd.read_csv(edgelist_path, sep='\t')

# Load ENSEMBL to gene symbol mapping file (ensure it has columns like 'ENSEMBL_ID' and 'Gene_Symbol')
mapping_df = pd.read_csv(mapping_path)

# Merge mapping with the edgelist based on gene1 and gene2
edgelist_df = edgelist_df.merge(mapping_df, how='left', left_on='gene1', right_on='ENSEMBL_ID')
edgelist_df = edgelist_df.merge(mapping_df, how='left', left_on='gene2', right_on='ENSEMBL_ID', suffixes=('_gene1', '_gene2'))

# Keep relevant columns
edgelist_df = edgelist_df[['gene1', 'gene2', 'Gene_Symbol_gene1', 'Gene_Symbol_gene2', 'module1', 'module2']]

# Filter rows where module1 is 'turquoise' and module2 is 'red'
filtered_df = edgelist_df[(edgelist_df['module1'] == 'turquoise') & (edgelist_df['module2'] == 'red')]

# Save the filtered edgelist with gene symbols
output_path = '/mnt/c/Users/agsko/dev/pcm/WGCNA/edgelist_filtered_turquoise_red.tsv'
filtered_df.to_csv(output_path, sep='\t', index=False)

print(f"Filtered edgelist with gene symbols saved to: {output_path}")
