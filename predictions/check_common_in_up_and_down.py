import pandas as pd

# Define the files and their source names
file_paths = {
    "/mnt/c/Users/agsko/dev/pcm/EdgeR/downregulated and upregulated genes/Ure_DMXAA_vs_PBS_RT_DMXAA_upregulated_genes.csv": "Ure_DMXAA_vs_PBS_RT_DMXAA",
    "/mnt/c/Users/agsko/dev/pcm/EdgeR/downregulated and upregulated genes/Ure_DMXAA_vs_Ure_RT_upregulated_genes.csv": "Ure_DMXAA_vs_Ure_RT",
    "/mnt/c/Users/agsko/dev/pcm/EdgeR/downregulated and upregulated genes/Ure_DMXAA_vs_Ure_upregulated_genes.csv": "Ure_DMXAA_vs_Ure",
    "/mnt/c/Users/agsko/dev/pcm/EdgeR/downregulated and upregulated genes/Ure_RT_DMXAA_vs_PBS_RT_DMXAA_upregulated_genes.csv": "Ure_RT_DMXAA_vs_PBS_RT_DMXAA",
    "/mnt/c/Users/agsko/dev/pcm/EdgeR/downregulated and upregulated genes/Ure_RT_DMXAA_vs_Ure_DMXAA_upregulated_genes.csv": "Ure_RT_DMXAA_vs_Ure_DMXAA",
    "/mnt/c/Users/agsko/dev/pcm/EdgeR/downregulated and upregulated genes/Ure_RT_DMXAA_vs_Ure_RT_upregulated_genes.csv": "Ure_RT_DMXAA_vs_Ure_RT",
    "/mnt/c/Users/agsko/dev/pcm/EdgeR/downregulated and upregulated genes/Ure_RT_DMXAA_vs_Ure_upregulated_genes.csv": "Ure_RT_DMXAA_vs_Ure",
    "/mnt/c/Users/agsko/dev/pcm/EdgeR/downregulated and upregulated genes/Ure_RT_vs_PBS_RT_DMXAA_upregulated_genes.csv": "Ure_RT_vs_PBS_RT_DMXAA",
    "/mnt/c/Users/agsko/dev/pcm/EdgeR/downregulated and upregulated genes/Ure_RT_vs_Ure_upregulated_genes.csv": "Ure_RT_vs_Ure",
    "/mnt/c/Users/agsko/dev/pcm/EdgeR/downregulated and upregulated genes/Ure_vs_PBS_RT_DMXAA_upregulated_genes.csv": "Ure_vs_PBS_RT_DMXAA"
}

# List to hold the dataframes
dfs = []

# Read each file and add a column indicating its source
for file_path, source in file_paths.items():
    df = pd.read_csv(file_path)
    df['Source'] = source  # Add a new column with the source name
    dfs.append(df)  # Append the dataframe to the list

# Concatenate all the dataframes
merged_df = pd.concat(dfs, ignore_index=True)

# Save the merged dataframe to a new file
output_file = "/mnt/c/Users/agsko/dev/pcm/predictions/merged_upregulated_genes.csv"
merged_df.to_csv(output_file, index=False)

print(f"Merged file saved to {output_file}")

import pandas as pd

# Define the file paths and source names for the downregulated genes
downregulated_file_paths = {
    "/mnt/c/Users/agsko/dev/pcm/EdgeR/downregulated and upregulated genes/Ure_DMXAA_vs_PBS_RT_DMXAA_downregulated_genes.csv": "Ure_DMXAA_vs_PBS_RT_DMXAA",
    "/mnt/c/Users/agsko/dev/pcm/EdgeR/downregulated and upregulated genes/Ure_DMXAA_vs_Ure_RT_downregulated_genes.csv": "Ure_DMXAA_vs_Ure_RT",
    "/mnt/c/Users/agsko/dev/pcm/EdgeR/downregulated and upregulated genes/Ure_DMXAA_vs_Ure_downregulated_genes.csv": "Ure_DMXAA_vs_Ure",
    "/mnt/c/Users/agsko/dev/pcm/EdgeR/downregulated and upregulated genes/Ure_RT_DMXAA_vs_PBS_RT_DMXAA_downregulated_genes.csv": "Ure_RT_DMXAA_vs_PBS_RT_DMXAA",
    "/mnt/c/Users/agsko/dev/pcm/EdgeR/downregulated and upregulated genes/Ure_RT_DMXAA_vs_Ure_DMXAA_downregulated_genes.csv": "Ure_RT_DMXAA_vs_Ure_DMXAA",
    "/mnt/c/Users/agsko/dev/pcm/EdgeR/downregulated and upregulated genes/Ure_RT_DMXAA_vs_Ure_RT_downregulated_genes.csv": "Ure_RT_DMXAA_vs_Ure_RT",
    "/mnt/c/Users/agsko/dev/pcm/EdgeR/downregulated and upregulated genes/Ure_RT_DMXAA_vs_Ure_downregulated_genes.csv": "Ure_RT_DMXAA_vs_Ure",
    "/mnt/c/Users/agsko/dev/pcm/EdgeR/downregulated and upregulated genes/Ure_RT_vs_PBS_RT_DMXAA_downregulated_genes.csv": "Ure_RT_vs_PBS_RT_DMXAA",
    "/mnt/c/Users/agsko/dev/pcm/EdgeR/downregulated and upregulated genes/Ure_RT_vs_Ure_downregulated_genes.csv": "Ure_RT_vs_Ure",
    "/mnt/c/Users/agsko/dev/pcm/EdgeR/downregulated and upregulated genes/Ure_vs_PBS_RT_DMXAA_downregulated_genes.csv": "Ure_vs_PBS_RT_DMXAA"
}

# List to hold the dataframes
downregulated_dfs = []

# Read each downregulated file and add a column indicating its source
for file_path, source in downregulated_file_paths.items():
    df = pd.read_csv(file_path)
    df['Source'] = source  # Add a new column with the source name
    downregulated_dfs.append(df)  # Append the dataframe to the list

# Concatenate all the dataframes
merged_downregulated_df = pd.concat(downregulated_dfs, ignore_index=True)

# Save the merged dataframe to a new file
downregulated_output_file = "/mnt/c/Users/agsko/dev/pcm/predictions/merged_downregulated_genes.csv"
merged_downregulated_df.to_csv(downregulated_output_file, index=False)

print(f"Downregulated merged file saved to {downregulated_output_file}")

import pandas as pd

# Define the file paths
common_peptides_file = "/mnt/c/Users/agsko/dev/pcm/predictions/common_peptides_in_different_samples.csv"
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
output_file = "/mnt/c/Users/agsko/dev/pcm/predictions/common_genes_in_peptides_upregulated_or_downregulated.csv"
common_genes.to_csv(output_file, index=False)

print(f"Common genes saved to {output_file}")

# Check if any genes from the common peptides file are present in the upregulated genes
common_genes_upregulated = pd.merge(common_peptides_df, upregulated_genes_df, on="Gene.name", how="inner")

# Save the result to a new file
output_file_upregulated = "/mnt/c/Users/agsko/dev/pcm/predictions/common_genes_in_peptides_upregulated_class1.csv"
common_genes_upregulated.to_csv(output_file_upregulated, index=False)

print(f"Common genes saved to {output_file_upregulated}")