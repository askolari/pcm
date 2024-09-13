import pandas as pd
import os

# Define the directory containing the merged files and where to save the output
merged_files_dir = '/mnt/c/Users/agsko/dev/pcm/merged_files/'
output_directory = '/mnt/c/Users/agsko/dev/pcm/WES_explor/common_gene_analysis'

# Load the final dataset into a pandas DataFrame
df = pd.read_csv(os.path.join(merged_files_dir, 'final_variants_mm10.csv'))

# Ensure the SYMBOL column is treated as a string and handle NaN values
df['SYMBOL'] = df['SYMBOL'].astype(str).fillna('')

# Extract unique genes (SYMBOL) for each group
ure_ir_genes = set(df[df['Sample_Group'] == 'Ure_IR']['SYMBOL'].unique())
ure_mock_ir_genes = set(df[df['Sample_Group'] == 'Ure_mock_IR']['SYMBOL'].unique())
ure_mock_ir_dmxaa_genes = set(df[df['Sample_Group'] == 'Ure_mock_IR_DMXAA']['SYMBOL'].unique())
pbs_ir_dmxaa_genes = set(df[df['Sample_Group'] == 'PBS_IR_DMXAA']['SYMBOL'].unique())
control_genes = set(df[df['Sample_Group'] == 'control']['SYMBOL'].unique())

# Remove any empty strings if present
ure_ir_genes = {gene for gene in ure_ir_genes if gene.strip()}
ure_mock_ir_genes = {gene for gene in ure_mock_ir_genes if gene.strip()}
ure_mock_ir_dmxaa_genes = {gene for gene in ure_mock_ir_dmxaa_genes if gene.strip()}
pbs_ir_dmxaa_genes = {gene for gene in pbs_ir_dmxaa_genes if gene.strip()}
control_genes = {gene for gene in control_genes if gene.strip()}

# 1. Common genes among all three Ure conditions
common_all_three_ure_conditions = ure_ir_genes.intersection(ure_mock_ir_genes).intersection(ure_mock_ir_dmxaa_genes)

# 2. Exclude genes that are present in PBS and control
exclusive_ure_genes = common_all_three_ure_conditions.difference(pbs_ir_dmxaa_genes).difference(control_genes)

# Save the exclusive genes to a CSV file
exclusive_genes_file = os.path.join(output_directory, 'exclusive_ure_genes_excluding_pbs_control.csv')
pd.DataFrame({'Exclusive_Genes': list(exclusive_ure_genes)}).to_csv(exclusive_genes_file, index=False)

print(f"Exclusive Ure condition genes (excluding PBS and control) saved to {exclusive_genes_file}")

import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

# Define the directories
merged_files_dir = '/mnt/c/Users/agsko/dev/pcm/merged_files/'
output_directory = '/mnt/c/Users/agsko/dev/pcm/WES_explor/common_gene_analysis/'

# Load the final dataset
df = pd.read_csv(os.path.join(merged_files_dir, 'final_variants_mm10.csv'))

# Load the exclusive genes list
exclusive_genes_file = os.path.join(output_directory, 'exclusive_ure_genes_excluding_pbs_control.csv')
exclusive_genes_df = pd.read_csv(exclusive_genes_file)
exclusive_genes = set(exclusive_genes_df['Exclusive_Genes'].dropna().unique())

# Filter the original dataframe to keep only the rows corresponding to the exclusive genes
filtered_df = df[df['SYMBOL'].isin(exclusive_genes)]

# Analyze the types of consequences associated with these genes

# Plot 1: Frequency of Consequence Types for Exclusive Genes
plt.figure(figsize=(12, 8))
sns.countplot(y='Consequence', data=filtered_df, order=filtered_df['Consequence'].value_counts().index, palette="coolwarm")
plt.title('Distribution of Consequences for Exclusive Ure Condition Genes')
plt.xlabel('Count')
plt.ylabel('Consequence')
plt.savefig(os.path.join(output_directory, 'exclusive_ure_genes_consequence_distribution.png'))
plt.show()

# Plot 2: Box plot for VAF across different Consequence Types
plt.figure(figsize=(12, 8))
sns.boxplot(x='Consequence', y='VAF', data=filtered_df)
plt.title('Variant Allele Frequencies by Consequence Type')
plt.xlabel('Consequence')
plt.ylabel('VAF')
plt.xticks(rotation=45)
plt.savefig(os.path.join(output_directory, 'exclusive_ure_genes_vaf_by_consequence.png'))
plt.show()

# Plot 3: Distribution of CSQ Impact with color differentiation
plt.figure(figsize=(10, 6))
sns.countplot(x='IMPACT', data=filtered_df, palette="Set2")
plt.title('Distribution of Impact for Exclusive Ure Condition Genes')
plt.xlabel('Impact')
plt.ylabel('Count')
plt.savefig(os.path.join(output_directory, 'exclusive_ure_genes_impact_distribution.png'))
plt.show()

# Save the filtered dataframe for further analysis
filtered_df.to_csv(os.path.join(output_directory, 'exclusive_ure_genes_filtered_data.csv'), index=False)

import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

# Define the directories
merged_files_dir = '/mnt/c/Users/agsko/dev/pcm/merged_files/'
output_directory = '/mnt/c/Users/agsko/dev/pcm/WES_explor/common_gene_analysis/'

# Load the final dataset
df = pd.read_csv(os.path.join(merged_files_dir, 'final_variants_mm10.csv'))

# Load the exclusive genes list
exclusive_genes_file = os.path.join(output_directory, 'exclusive_ure_genes_excluding_pbs_control.csv')
exclusive_genes_df = pd.read_csv(exclusive_genes_file)
exclusive_genes = set(exclusive_genes_df['Exclusive_Genes'].dropna().unique())

# Filter the original dataframe to keep only the rows corresponding to the exclusive genes
filtered_df = df[df['SYMBOL'].isin(exclusive_genes)]

# Analyze the types of consequences associated with these genes

# Plot 1: Frequency of Consequence Types for Exclusive Genes
plt.figure(figsize=(12, 8))
sns.countplot(y='Consequence', data=filtered_df, order=filtered_df['Consequence'].value_counts().index, palette="coolwarm")
plt.title('Distribution of Consequences for Exclusive Ure Condition Genes')
plt.xlabel('Count')
plt.ylabel('Consequence')
plt.savefig(os.path.join(output_directory, 'exclusive_ure_genes_consequence_distribution.png'))
plt.show()

# Plot 2: Box plot for VAF across different Consequence Types
plt.figure(figsize=(12, 8))
sns.boxplot(x='Consequence', y='VAF', data=filtered_df)
plt.title('Variant Allele Frequencies by Consequence Type')
plt.xlabel('Consequence')
plt.ylabel('VAF')
plt.xticks(rotation=45)
plt.savefig(os.path.join(output_directory, 'exclusive_ure_genes_vaf_by_consequence.png'))
plt.show()

# Plot 3: Distribution of CSQ Impact with color differentiation
plt.figure(figsize=(10, 6))
sns.countplot(x='IMPACT', data=filtered_df, palette="Set2")
plt.title('Distribution of Impact for Exclusive Ure Condition Genes')
plt.xlabel('Impact')
plt.ylabel('Count')
plt.savefig(os.path.join(output_directory, 'exclusive_ure_genes_impact_distribution.png'))
plt.show()

# Save the filtered dataframe for further analysis
filtered_df.to_csv(os.path.join(output_directory, 'exclusive_ure_genes_filtered_data.csv'), index=False)
