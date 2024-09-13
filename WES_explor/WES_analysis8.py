import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

# Define directories
merged_files_dir = '/mnt/c/Users/agsko/dev/pcm/merged_files/'
output_directory = '/mnt/c/Users/agsko/dev/pcm/WES_explor/common_gene_analysis'
individual_gene_dir = os.path.join(output_directory, 'individual_gene_analysis')
dataset_output_path = os.path.join(output_directory, 'exclusive_gene_analysis_dataset.csv')

# Ensure the output directory exists
if not os.path.exists(individual_gene_dir):
    os.makedirs(individual_gene_dir)

# Load the final dataset
df = pd.read_csv(os.path.join(merged_files_dir, 'final_variants_mm10.csv'))

# Load the exclusive genes list
exclusive_genes_file = os.path.join(output_directory, 'exclusive_ure_genes_excluding_pbs_control.csv')
exclusive_genes_df = pd.read_csv(exclusive_genes_file)
exclusive_genes = set(exclusive_genes_df['Exclusive_Genes'].dropna().unique())

# Filter the original dataframe to keep only the rows corresponding to the exclusive genes
filtered_df = df[df['SYMBOL'].isin(exclusive_genes)]

# Create a list to hold all relevant information for the new dataset
all_gene_data = []

# Iterate through each exclusive gene and generate plots
for gene in exclusive_genes:
    gene_df = filtered_df[filtered_df['SYMBOL'] == gene]
    
    # Collect relevant information for the dataset
    for _, row in gene_df.iterrows():
        all_gene_data.append({
            'SYMBOL': gene,
            'Sample_Group': row['Sample_Group'],
            'VAF': row['VAF'],
            'IMPACT': row['IMPACT'],
            'Consequence': row['Consequence']
        })
    
    # Plot 1: VAF distribution by Sample Group for this gene
    plt.figure(figsize=(12, 8))
    sns.boxplot(x='Sample_Group', y='VAF', data=gene_df, palette="Set2")
    plt.title(f'VAF Distribution for {gene} by Sample Group')
    plt.xlabel('Sample Group')
    plt.ylabel('VAF')
    plt.xticks(rotation=45)
    plt.show()  # Display the plot
    plt.savefig(os.path.join(individual_gene_dir, f'{gene}_vaf_distribution_by_sample_group.png'))
    plt.close()
    
    # Plot 2: Impact distribution for this gene across sample groups
    plt.figure(figsize=(12, 8))
    sns.countplot(x='Sample_Group', hue='IMPACT', data=gene_df, palette="Set2")
    plt.title(f'Impact Distribution for {gene} by Sample Group')
    plt.xlabel('Sample Group')
    plt.ylabel('Count')
    plt.xticks(rotation=45)
    plt.show()  # Display the plot
    plt.savefig(os.path.join(individual_gene_dir, f'{gene}_impact_distribution_by_sample_group.png'))
    plt.close()
    
    # Plot 3: Consequence distribution for this gene across sample groups
    plt.figure(figsize=(12, 8))
    sns.countplot(x='Sample_Group', hue='Consequence', data=gene_df, palette="coolwarm")
    plt.title(f'Consequence Distribution for {gene} by Sample Group')
    plt.xlabel('Sample Group')
    plt.ylabel('Count')
    plt.xticks(rotation=45)
    plt.show()  # Display the plot
    plt.savefig(os.path.join(individual_gene_dir, f'{gene}_consequence_distribution_by_sample_group.png'))
    plt.close()

    print(f'Plots for {gene} saved successfully in {individual_gene_dir}')

# Convert the list of all gene data to a DataFrame
all_gene_data_df = pd.DataFrame(all_gene_data)

# Save the new dataset to a CSV file
all_gene_data_df.to_csv(dataset_output_path, index=False)

print(f'All gene-specific analyses have been saved in {individual_gene_dir}')
print(f'New dataset with all relevant information saved to {dataset_output_path}')
