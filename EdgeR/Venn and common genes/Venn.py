import pandas as pd
from matplotlib_venn import venn3
import matplotlib.pyplot as plt

# Adjust the path to the Linux-compatible format for downregulated gene sets
ure_rt_down_path = "/mnt/c/Users/agsko/dev/pcm/EdgeR/Ure_RT_vs_PBS_RT_DMXAA_downregulated_genes.csv"
ure_dmxaa_down_path = "/mnt/c/Users/agsko/dev/pcm/EdgeR/Ure_DMXAA_vs_PBS_RT_DMXAA_downregulated_genes.csv"
ure_rt_dmxaa_down_path = "/mnt/c/Users/agsko/dev/pcm/EdgeR/Ure_RT_DMXAA_vs_PBS_RT_DMXAA_downregulated_genes.csv"

# Load downregulated gene sets (assuming the gene ID is in 'ID' and gene name in 'Gene.name')
ure_rt_down = pd.read_csv(ure_rt_down_path).dropna(subset=['ID', 'Gene.name'])
ure_dmxaa_down = pd.read_csv(ure_dmxaa_down_path).dropna(subset=['ID', 'Gene.name'])
ure_rt_dmxaa_down = pd.read_csv(ure_rt_dmxaa_down_path).dropna(subset=['ID', 'Gene.name'])

# Create a Venn diagram for downregulated genes (using 'ID' column for comparison)
plt.figure(figsize=(8, 8))
venn3([set(ure_rt_down['ID']), set(ure_dmxaa_down['ID']), set(ure_rt_dmxaa_down['ID'])], 
      set_labels=('Ure_RT', 'Ure_DMXAA', 'Ure_RT_DMXAA'))

# Add a title
plt.title("Venn Diagram of Downregulated Genes")
plt.show()

# Find common downregulated genes between all three conditions by gene ID
common_all_three_down_ids = set(ure_rt_down['ID']).intersection(set(ure_dmxaa_down['ID']), set(ure_rt_dmxaa_down['ID']))
common_all_three_down = ure_rt_down[ure_rt_down['ID'].isin(common_all_three_down_ids)]

# Find common downregulated genes between each pair of conditions by gene ID
common_ure_rt_ure_dmxaa_down_ids = set(ure_rt_down['ID']).intersection(set(ure_dmxaa_down['ID']))
common_ure_rt_ure_dmxaa_down = ure_rt_down[ure_rt_down['ID'].isin(common_ure_rt_ure_dmxaa_down_ids)]

common_ure_dmxaa_ure_rt_dmxaa_down_ids = set(ure_dmxaa_down['ID']).intersection(set(ure_rt_dmxaa_down['ID']))
common_ure_dmxaa_ure_rt_dmxaa_down = ure_dmxaa_down[ure_dmxaa_down['ID'].isin(common_ure_dmxaa_ure_rt_dmxaa_down_ids)]

common_ure_rt_ure_rt_dmxaa_down_ids = set(ure_rt_down['ID']).intersection(set(ure_rt_dmxaa_down['ID']))
common_ure_rt_ure_rt_dmxaa_down = ure_rt_down[ure_rt_down['ID'].isin(common_ure_rt_ure_rt_dmxaa_down_ids)]

# Find unique downregulated genes in each condition by gene ID
unique_ure_rt_down_ids = set(ure_rt_down['ID']) - (set(ure_dmxaa_down['ID']).union(set(ure_rt_dmxaa_down['ID'])))
unique_ure_rt_down = ure_rt_down[ure_rt_down['ID'].isin(unique_ure_rt_down_ids)]

unique_ure_dmxaa_down_ids = set(ure_dmxaa_down['ID']) - (set(ure_rt_down['ID']).union(set(ure_rt_dmxaa_down['ID'])))
unique_ure_dmxaa_down = ure_dmxaa_down[ure_dmxaa_down['ID'].isin(unique_ure_dmxaa_down_ids)]

unique_ure_rt_dmxaa_down_ids = set(ure_rt_dmxaa_down['ID']) - (set(ure_rt_down['ID']).union(set(ure_dmxaa_down['ID'])))
unique_ure_rt_dmxaa_down = ure_rt_dmxaa_down[ure_rt_dmxaa_down['ID'].isin(unique_ure_rt_dmxaa_down_ids)]

# Save the common and unique sets (with both gene IDs and gene names) to CSV files
common_all_three_down.to_csv("/mnt/c/Users/agsko/dev/pcm/EdgeR/common_all_ure_RT_ure_DMXAA_ure_RT_DMXAA_downregulated_with_names.csv", index=False)
common_ure_rt_ure_dmxaa_down.to_csv("/mnt/c/Users/agsko/dev/pcm/EdgeR/common_ure_rt_ure_dmxaa_genes_downregulated_with_names.csv", index=False)
common_ure_dmxaa_ure_rt_dmxaa_down.to_csv("/mnt/c/Users/agsko/dev/pcm/EdgeR/common_ure_dmxaa_ure_rt_dmxaa_genes_downregulated_with_names.csv", index=False)
common_ure_rt_ure_rt_dmxaa_down.to_csv("/mnt/c/Users/agsko/dev/pcm/EdgeR/common_ure_rt_ure_rt_dmxaa_genes_downregulated_with_names.csv", index=False)

unique_ure_rt_down.to_csv("/mnt/c/Users/agsko/dev/pcm/EdgeR/unique_ure_rt_genes_downregulated_with_names.csv", index=False)
unique_ure_dmxaa_down.to_csv("/mnt/c/Users/agsko/dev/pcm/EdgeR/unique_ure_dmxaa_genes_downregulated_with_names.csv", index=False)
unique_ure_rt_dmxaa_down.to_csv("/mnt/c/Users/agsko/dev/pcm/EdgeR/unique_ure_rt_dmxaa_genes_downregulated_with_names.csv", index=False)


# Repeat for upregulated genes
ure_rt_up_path = "/mnt/c/Users/agsko/dev/pcm/EdgeR/Ure_RT_vs_PBS_RT_DMXAA_upregulated_genes.csv"
ure_dmxaa_up_path = "/mnt/c/Users/agsko/dev/pcm/EdgeR/Ure_DMXAA_vs_PBS_RT_DMXAA_upregulated_genes.csv"
ure_rt_dmxaa_up_path = "/mnt/c/Users/agsko/dev/pcm/EdgeR/Ure_RT_DMXAA_vs_PBS_RT_DMXAA_upregulated_genes.csv"

# Load upregulated gene sets (assuming the gene ID is in 'ID' and gene name in 'Gene.name')
ure_rt_up = pd.read_csv(ure_rt_up_path).dropna(subset=['ID', 'Gene.name'])
ure_dmxaa_up = pd.read_csv(ure_dmxaa_up_path).dropna(subset=['ID', 'Gene.name'])
ure_rt_dmxaa_up = pd.read_csv(ure_rt_dmxaa_up_path).dropna(subset=['ID', 'Gene.name'])

# Create a Venn diagram for upregulated genes (using 'ID' column for comparison)
plt.figure(figsize=(8, 8))
venn3([set(ure_rt_up['ID']), set(ure_dmxaa_up['ID']), set(ure_rt_dmxaa_up['ID'])], 
      set_labels=('Ure_RT', 'Ure_DMXAA', 'Ure_RT_DMXAA'))

# Add a title
plt.title("Venn Diagram of Upregulated Genes")
plt.show()

# Find common upregulated genes between all three conditions by gene ID
common_all_three_up_ids = set(ure_rt_up['ID']).intersection(set(ure_dmxaa_up['ID']), set(ure_rt_dmxaa_up['ID']))
common_all_three_up = ure_rt_up[ure_rt_up['ID'].isin(common_all_three_up_ids)]

# Find common upregulated genes between each pair of conditions by gene ID
common_ure_rt_ure_dmxaa_up_ids = set(ure_rt_up['ID']).intersection(set(ure_dmxaa_up['ID']))
common_ure_rt_ure_dmxaa_up = ure_rt_up[ure_rt_up['ID'].isin(common_ure_rt_ure_dmxaa_up_ids)]

common_ure_dmxaa_ure_rt_dmxaa_up_ids = set(ure_dmxaa_up['ID']).intersection(set(ure_rt_dmxaa_up['ID']))
common_ure_dmxaa_ure_rt_dmxaa_up = ure_dmxaa_up[ure_dmxaa_up['ID'].isin(common_ure_dmxaa_ure_rt_dmxaa_up_ids)]

common_ure_rt_ure_rt_dmxaa_up_ids = set(ure_rt_up['ID']).intersection(set(ure_rt_dmxaa_up['ID']))
common_ure_rt_ure_rt_dmxaa_up = ure_rt_up[ure_rt_up['ID'].isin(common_ure_rt_ure_rt_dmxaa_up_ids)]

# Find unique upregulated genes in each condition by gene ID
unique_ure_rt_up_ids = set(ure_rt_up['ID']) - (set(ure_dmxaa_up['ID']).union(set(ure_rt_dmxaa_up['ID'])))
unique_ure_rt_up = ure_rt_up[ure_rt_up['ID'].isin(unique_ure_rt_up_ids)]

unique_ure_dmxaa_up_ids = set(ure_dmxaa_up['ID']) - (set(ure_rt_up['ID']).union(set(ure_rt_dmxaa_up['ID'])))
unique_ure_dmxaa_up = ure_dmxaa_up[ure_dmxaa_up['ID'].isin(unique_ure_dmxaa_up_ids)]

unique_ure_rt_dmxaa_up_ids = set(ure_rt_dmxaa_up['ID']) - (set(ure_rt_up['ID']).union(set(ure_dmxaa_up['ID'])))
unique_ure_rt_dmxaa_up = ure_rt_dmxaa_up[ure_rt_dmxaa_up['ID'].isin(unique_ure_rt_dmxaa_up_ids)]

# Save the common and unique sets (with both gene IDs and gene names) to CSV files
common_all_three_up.to_csv("/mnt/c/Users/agsko/dev/pcm/EdgeR/common_all_ure_RT_ure_DMXAA_ure_RT_DMXAA_upregulated_with_names.csv", index=False)
common_ure_rt_ure_dmxaa_up.to_csv("/mnt/c/Users/agsko/dev/pcm/EdgeR/common_ure_rt_ure_dmxaa_genes_upregulated_with_names.csv", index=False)
common_ure_dmxaa_ure_rt_dmxaa_up.to_csv("/mnt/c/Users/agsko/dev/pcm/EdgeR/common_ure_dmxaa_ure_rt_dmxaa_genes_upregulated_with_names.csv", index=False)
common_ure_rt_ure_rt_dmxaa_up.to_csv("/mnt/c/Users/agsko/dev/pcm/EdgeR/common_ure_rt_ure_rt_dmxaa_genes_upregulated_with_names.csv", index=False)

unique_ure_rt_up.to_csv("/mnt/c/Users/agsko/dev/pcm/EdgeR/unique_ure_rt_genes_upregulated_with_names.csv", index=False)
unique_ure_dmxaa_up.to_csv("/mnt/c/Users/agsko/dev/pcm/EdgeR/unique_ure_dmxaa_genes_upregulated_with_names.csv", index=False)
unique_ure_rt_dmxaa_up.to_csv("/mnt/c/Users/agsko/dev/pcm/EdgeR/unique_ure_rt_dmxaa_genes_upregulated_with_names.csv", index=False)
