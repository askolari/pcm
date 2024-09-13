import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Define the directory containing the final dataset and output directory for images
data_dir = '/mnt/c/Users/agsko/dev/pcm/merged_files/'
output_dir = '/mnt/c/Users/agsko/dev/pcm/WES_explor/top_genes_analysis'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Load the final dataset
df = pd.read_csv(os.path.join(data_dir, 'final_variants_mm10.csv'))

# Identify the top genes by frequency of variants in each group
top_n = 10  # Adjust as needed
grouped = df.groupby(['Sample_Group', 'SYMBOL']).size().reset_index(name='counts')
top_genes_by_group = grouped.groupby('Sample_Group').apply(lambda x: x.nlargest(top_n, 'counts')).reset_index(drop=True)

# Save the top genes by group
top_genes_file = os.path.join(output_dir, 'top_genes_by_group.csv')
top_genes_by_group.to_csv(top_genes_file, index=False)
print(f"Top genes by group saved to {top_genes_file}")

# Identify common genes among the groups
ure_ir_genes = set(df[df['Sample_Group'] == 'Ure_IR']['SYMBOL'].unique())
ure_mock_ir_genes = set(df[df['Sample_Group'] == 'Ure_mock_IR']['SYMBOL'].unique())
ure_mock_ir_dmxaa_genes = set(df[df['Sample_Group'] == 'Ure_mock_IR_DMXAA']['SYMBOL'].unique())

# Common genes in all Ure conditions
common_genes_all_ure = ure_ir_genes & ure_mock_ir_genes & ure_mock_ir_dmxaa_genes

# Save common genes
common_genes_file = os.path.join(output_dir, 'common_genes_all_ure.csv')
pd.DataFrame(list(common_genes_all_ure), columns=['SYMBOL']).to_csv(common_genes_file, index=False)
print(f"Common genes among all Ure conditions saved to {common_genes_file}")

# Further analysis and visualization on top genes and common genes
def visualize_genes(genes, df, prefix):
    # Filter the dataframe for the selected genes
    df_filtered = df[df['SYMBOL'].isin(genes)]

    # Plot: Consequence distribution
    plt.figure(figsize=(14, 8))
    sns.countplot(y='Consequence', data=df_filtered, order=df_filtered['Consequence'].value_counts().index, palette="viridis")
    plt.title(f'Distribution of Consequences for {prefix}')
    plt.xlabel('Count')
    plt.ylabel('Consequence')
    plt.savefig(os.path.join(output_dir, f'{prefix}_consequence_distribution.png'))
    plt.show()

    # Plot: Impact distribution
    plt.figure(figsize=(10, 6))
    sns.countplot(x='IMPACT', data=df_filtered, palette="Set2")
    plt.title(f'Distribution of IMPACT for {prefix}')
    plt.xlabel('Impact')
    plt.ylabel('Count')
    plt.savefig(os.path.join(output_dir, f'{prefix}_impact_distribution.png'))
    plt.show()

    # Plot: VAF distribution by Consequence
    plt.figure(figsize=(14, 8))
    sns.boxplot(x='VAF', y='Consequence', data=df_filtered, palette="coolwarm")
    plt.title(f'VAF Distribution by Consequence for {prefix}')
    plt.xlabel('Variant Allele Frequency (VAF)')
    plt.ylabel('Consequence')
    plt.savefig(os.path.join(output_dir, f'{prefix}_vaf_distribution_by_consequence.png'))
    plt.show()

# Visualize top genes by group
for group in top_genes_by_group['Sample_Group'].unique():
    top_genes_in_group = top_genes_by_group[top_genes_by_group['Sample_Group'] == group]['SYMBOL']
    visualize_genes(top_genes_in_group, df, f'top_genes_{group}')

# Visualize common genes among all Ure conditions
visualize_genes(common_genes_all_ure, df, 'common_genes_all_ure')
