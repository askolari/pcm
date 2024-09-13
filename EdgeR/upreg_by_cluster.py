import pandas as pd
from matplotlib_venn import venn2
import matplotlib.pyplot as plt

# Load the two CSV files containing upregulated genes
path1 = "/mnt/c/Users/agsko/dev/pcm/EdgeR/downregulated and upregulated genes/Ure_RT_vs_PBS_RT_DMXAA_upregulated_genes.csv"
path2 = "/mnt/c/Users/agsko/dev/pcm/EdgeR/dge_cluster4_pbs_vs_ure_rt_upregulated_with_entrez.csv"

# Load the data from the CSV files
geneset1 = pd.read_csv(path1)
geneset2 = pd.read_csv(path2)

# Extract the 'Gene.name' columns from both datasets
geneset1_genes = set(geneset1['Gene.name'])
geneset2_genes = set(geneset2['Gene.name'])

# Create a Venn diagram comparing the two sets of genes
plt.figure(figsize=(8, 8))
venn2([geneset1_genes, geneset2_genes], set_labels=('Ure_RT_vs_PBS_RT_DMXAA', 'cluster4_pbs_vs_ure_rt'))
# Display the Venn diagram
plt.title('Venn Diagram of Upregulated Genes')
plt.show()

import os

# Extract the file names (without extensions) from the paths for use in the output file name
contrast1 = os.path.splitext(os.path.basename(path1))[0]
contrast2 = os.path.splitext(os.path.basename(path2))[0]

# Load the data from the CSV files
geneset1 = pd.read_csv(path1)
geneset2 = pd.read_csv(path2)

# Merge the datasets based on 'Gene.name' using an outer join to include all genes
merged_genes = pd.merge(geneset1, geneset2, on='Gene.name', how='outer', suffixes=('_Ure_RT_vs_PBS', '_Cluster1'))

# Create a dynamic file name based on the contrasts
output_filename = f"/mnt/c/Users/agsko/dev/pcm/EdgeR/merged_{contrast1}_and_{contrast2}_upregulated_genes.csv"

# Save the merged dataset to a new CSV file with the dynamic name
merged_genes.to_csv(output_filename, index=False)
