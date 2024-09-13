import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2, venn3
import os

# Define the directory containing the merged files and where to save the output
merged_files_dir = '/mnt/c/Users/agsko/dev/pcm/merged_files/'
output_directory = '/mnt/c/Users/agsko/dev/pcm/WES_explor/additional_comparisons'

# Ensure the output directory exists
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Load the final dataset into a pandas DataFrame
df = pd.read_csv(os.path.join(merged_files_dir, 'final_variants_mm10.csv'))

# Filter for high and moderate impact variants
significant_impact_df = df[df['IMPACT'].isin(['HIGH', 'MODERATE'])]

# 1. Comparison between Ure_IR, Ure_mock_IR, and Ure_mock_IR_DMXAA
comparison_groups = ['Ure_IR', 'Ure_mock_IR', 'Ure_mock_IR_DMXAA']
comparison_df = significant_impact_df[significant_impact_df['Sample_Group'].isin(comparison_groups)]

# Save the filtered dataframe for Ure comparisons
comparison_file_path = os.path.join(merged_files_dir, 'Ure_comparison_significant_variants.csv')
comparison_df.to_csv(comparison_file_path, index=False)
print(f"Significant variants for Ure_IR, Ure_mock_IR, and Ure_mock_IR_DMXAA saved to {comparison_file_path}")

# Compare the distribution of high and moderate impact variants across these groups
plt.figure(figsize=(12, 6))
sns.countplot(x='IMPACT', hue='Sample_Group', data=comparison_df, palette="Set2")
plt.title('Significant Variant Impact Distribution between Ure_IR, Ure_mock_IR, and Ure_mock_IR_DMXAA')
plt.xlabel('Impact')
plt.ylabel('Count')
plt.savefig(os.path.join(output_directory, 'significant_impact_ure_comparison.png'))
plt.show()

# Compare specific consequences between Ure_IR, Ure_mock_IR, and Ure_mock_IR_DMXAA
consequences_comparison_df = pd.crosstab(comparison_df['Consequence'], comparison_df['Sample_Group'])

plt.figure(figsize=(14, 10))
sns.heatmap(consequences_comparison_df, annot=True, cmap="YlGnBu", fmt='d')
plt.title('Consequence Comparison between Ure_IR, Ure_mock_IR, and Ure_mock_IR_DMXAA (Significant Variants)')
plt.savefig(os.path.join(output_directory, 'consequence_comparison_ure_conditions.png'))
plt.show()

# Compare VAF distribution between these groups for significant variants
plt.figure(figsize=(10, 6))
sns.boxplot(x='Sample_Group', y='VAF', data=comparison_df, palette="Set3")
plt.title('VAF Distribution between Ure_IR, Ure_mock_IR, and Ure_mock_IR_DMXAA (Significant Variants)')
plt.savefig(os.path.join(output_directory, 'vaf_distribution_ure_conditions.png'))
plt.show()

# 2. Analysis of PBS_IR_DMXAA group
pbs_group_df = significant_impact_df[significant_impact_df['Sample_Group'] == 'PBS_IR_DMXAA']

# Save the PBS_IR_DMXAA significant variants to a new file
pbs_file_path = os.path.join(merged_files_dir, 'PBS_IR_DMXAA_significant_variants.csv')
pbs_group_df.to_csv(pbs_file_path, index=False)
print(f"Significant variants for PBS_IR_DMXAA saved to {pbs_file_path}")

# Distribution of Variant Types for PBS_IR_DMXAA
sns.countplot(x='Variant_Type', data=pbs_group_df)
plt.title('Distribution of Variant Types for PBS_IR_DMXAA (Significant Variants)')
plt.savefig(os.path.join(output_directory, 'pbs_impact_variant_type_distribution.png'))
plt.show()

# Consequence distribution for PBS_IR_DMXAA
sns.countplot(y='Consequence', data=pbs_group_df, order=pbs_group_df['Consequence'].value_counts().index)
plt.title('Distribution of Consequences for PBS_IR_DMXAA (Significant Variants)')
plt.savefig(os.path.join(output_directory, 'pbs_consequence_distribution.png'))
plt.show()

# Top genes in PBS_IR_DMXAA
top_genes_pbs = pbs_group_df['Gene'].value_counts().head(10)
plt.figure(figsize=(12, 6))
sns.barplot(x=top_genes_pbs.values, y=top_genes_pbs.index, palette="Blues_d")
plt.title('Top 10 Genes in PBS_IR_DMXAA (Significant Variants)')
plt.xlabel('Count')
plt.ylabel('Gene')
plt.savefig(os.path.join(output_directory, 'top_10_genes_pbs_ir_dmxaa.png'))
plt.show()

# Venn Diagram: PBS_IR_DMXAA vs Ure conditions
ure_conditions = set(df[df['Sample_Group'].isin(comparison_groups)]['SYMBOL'].unique())
pbs_genes = set(pbs_group_df['SYMBOL'].unique())

plt.figure(figsize=(8, 8))
venn2([ure_conditions, pbs_genes],
      set_labels=('Ure Conditions', 'PBS_IR_DMXAA'),
      set_colors=('skyblue', 'lightgreen'))
plt.title("Venn Diagram: Ure Conditions vs PBS_IR_DMXAA")
plt.savefig(os.path.join(output_directory, 'venn_ure_vs_pbs.png'))
plt.show()

# Unique significant variants in PBS_IR_DMXAA
unique_pbs_genes = pbs_genes.difference(ure_conditions)

print(f"Number of unique genes in PBS_IR_DMXAA: {len(unique_pbs_genes)}")
print("Unique genes in PBS_IR_DMXAA:")
for gene in unique_pbs_genes:
    print(gene)

# Save unique PBS_IR_DMXAA genes to a file
unique_pbs_file_path = os.path.join(output_directory, 'unique_pbs_ir_dmxaa_genes.csv')
pd.DataFrame({'Unique_Genes': list(unique_pbs_genes)}).to_csv(unique_pbs_file_path, index=False)
print(f"Unique PBS_IR_DMXAA genes saved to {unique_pbs_file_path}")
