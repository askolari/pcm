# Import necessary libraries
import pandas as pd
from matplotlib_venn import venn2
import matplotlib.pyplot as plt

# File paths
file1 = '/mnt/c/Users/agsko/dev/pcm/EdgeR/dge_cluster3_upregulated.csv'
file2 = '/mnt/c/Users/agsko/dev/pcm/EdgeR/downregulated and upregulated genes/Ure_DMXAA_vs_Ure_upregulated_genes.csv'

# Load the data
data1 = pd.read_csv(file1)
data2 = pd.read_csv(file2)

# Assuming the gene column is named 'Gene', modify if necessary
genes1 = set(data1['Gene'])
genes2 = set(data2['Gene'])

# Create the Venn diagram
plt.figure(figsize=(8, 8))
venn2([genes1, genes2], set_labels=('dge_cluster3_upregulated', 'Ure_DMXAA_vs_Ure_upregulated'))

# Show the plot
plt.title('Venn Diagram of Upregulated Genes')
plt.show()
