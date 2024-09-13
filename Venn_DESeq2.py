import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

# Define the directory where your files are located
directory = "/mnt/c/Users/agsko/dev/pcm/DESeq2"

# Define the file paths for all the CSV files
file_paths = {
    "Ure_vs_Ure_DMXAA": os.path.join(directory, "significant_genes_Ure_vs_Ure_DMXAA.csv"),
    "Ure_vs_Ure_RT": os.path.join(directory, "significant_genes_Ure_vs_Ure_RT.csv"),
    "Ure_vs_Ure_RT_DMXAA": os.path.join(directory, "significant_genes_Ure_vs_Ure_RT_DMXAA.csv"),
    "Ure_DMXAA_vs_PBS_RT_DMXAA": os.path.join(directory, "significant_genes_Ure_DMXAA_vs_PBS_RT_DMXAA.csv"),
    "Ure_DMXAA_vs_Ure_RT": os.path.join(directory, "significant_genes_Ure_DMXAA_vs_Ure_RT.csv"),
    "Ure_DMXAA_vs_Ure_RT_DMXAA": os.path.join(directory, "significant_genes_Ure_DMXAA_vs_Ure_RT_DMXAA.csv"),
    "Ure_RT_DMXAA_vs_PBS_RT_DMXAA": os.path.join(directory, "significant_genes_Ure_RT_DMXAA_vs_PBS_RT_DMXAA.csv"),
    "Ure_RT_vs_PBS_RT_DMXAA": os.path.join(directory, "significant_genes_Ure_RT_vs_PBS_RT_DMXAA.csv"),
    "Ure_RT_vs_Ure_RT_DMXAA": os.path.join(directory, "significant_genes_Ure_RT_vs_Ure_RT_DMXAA.csv"),
    "Ure_vs_PBS_RT_DMXAA": os.path.join(directory, "significant_genes_Ure_vs_PBS_RT_DMXAA.csv")
}

# Create a dictionary to store the gene lists from each file
gene_lists = {}

# Load each file and extract the first unnamed column (gene identifiers)
for contrast, file_path in file_paths.items():
    try:
        data = pd.read_csv(file_path)
        gene_lists[contrast] = data.iloc[:, 0].dropna().tolist()  # Use the first column for gene IDs
        print(f"Loaded {contrast} with {len(gene_lists[contrast])} genes.")
    except Exception as e:
        print(f"Error loading {contrast}: {e}")

# Convert the relevant lists to sets for the Venn diagram
set_Ure_vs_Ure_RT = set(gene_lists['Ure_vs_Ure_RT'])
set_Ure_RT_vs_Ure_RT_DMXAA = set(gene_lists['Ure_RT_vs_Ure_RT_DMXAA'])
set_Ure_RT_DMXAA_vs_PBS_RT_DMXAA = set(gene_lists['Ure_RT_DMXAA_vs_PBS_RT_DMXAA'])

# Create the Venn diagram
plt.figure(figsize=(5, 8))
venn = venn3([set_Ure_vs_Ure_RT, set_Ure_RT_vs_Ure_RT_DMXAA, set_Ure_RT_DMXAA_vs_PBS_RT_DMXAA],
             set_labels=('Ure vs Ure_RT', 'Ure_RT vs Ure_RT_DMXAA', 'Ure_RT_DMXAA vs PBS_RT_DMXAA'))

# Display the plot
plt.title("Venn Diagram of Gene Sets")
plt.show()

# Save the plot to a file
plt.savefig("venn_diagram_ure_conditions.png")

# Export common and unique genes
common_genes = set_Ure_vs_Ure_RT & set_Ure_RT_vs_Ure_RT_DMXAA & set_Ure_RT_DMXAA_vs_PBS_RT_DMXAA
unique_Ure_vs_Ure_RT = set_Ure_vs_Ure_RT - (set_Ure_RT_vs_Ure_RT_DMXAA | set_Ure_RT_DMXAA_vs_PBS_RT_DMXAA)
unique_Ure_RT_vs_Ure_RT_DMXAA = set_Ure_RT_vs_Ure_RT_DMXAA - (set_Ure_vs_Ure_RT | set_Ure_RT_DMXAA_vs_PBS_RT_DMXAA)
unique_Ure_RT_DMXAA_vs_PBS_RT_DMXAA = set_Ure_RT_DMXAA_vs_PBS_RT_DMXAA - (set_Ure_vs_Ure_RT | set_Ure_RT_vs_Ure_RT_DMXAA)

# Create DataFrames for export
common_genes_df = pd.DataFrame(list(common_genes), columns=["Gene"])
unique_Ure_vs_Ure_RT_df = pd.DataFrame(list(unique_Ure_vs_Ure_RT), columns=["Gene"])
unique_Ure_RT_vs_Ure_RT_DMXAA_df = pd.DataFrame(list(unique_Ure_RT_vs_Ure_RT_DMXAA), columns=["Gene"])
unique_Ure_RT_DMXAA_vs_PBS_RT_DMXAA_df = pd.DataFrame(list(unique_Ure_RT_DMXAA_vs_PBS_RT_DMXAA), columns=["Gene"])

# Export to CSV files
common_genes_df.to_csv("common_genes.csv", index=False)
unique_Ure_vs_Ure_RT_df.to_csv("unique_genes_Ure_vs_Ure_RT.csv", index=False)
unique_Ure_RT_vs_Ure_RT_DMXAA_df.to_csv("unique_genes_Ure_RT_vs_Ure_RT_DMXAA.csv", index=False)
unique_Ure_RT_DMXAA_vs_PBS_RT_DMXAA_df.to_csv("unique_genes_Ure_RT_DMXAA_vs_PBS_RT_DMXAA.csv", index=False)

print("Common and unique genes exported successfully!")