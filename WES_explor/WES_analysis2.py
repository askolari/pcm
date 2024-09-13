import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2, venn3
import os

# Specify the directory where you want to save the images
output_directory = '/mnt/c/Users/agsko/dev/pcm/WES_explor'

# Load your CSV file into a pandas DataFrame
df = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/merged_files/final_variants_mm10.csv')

# Convert the 'SYMBOL' arrays to sets
ure_ir_genes = set(df[df['Sample_Group'] == 'Ure_IR']['SYMBOL'].unique())
ure_mock_ir_genes = set(df[df['Sample_Group'] == 'Ure_mock_IR']['SYMBOL'].unique())
ure_mock_ir_dmxaa_genes = set(df[df['Sample_Group'] == 'Ure_mock_IR_DMXAA']['SYMBOL'].unique())
pbs_ir_dmxaa_genes = set(df[df['Sample_Group'] == 'PBS_IR_DMXAA']['SYMBOL'].unique())
control_genes = set(df[df['Sample_Group'] == 'control']['SYMBOL'].unique())

# Find common genes in Ure_IR and Ure_mock_IR
common_ure_genes = ure_ir_genes.intersection(ure_mock_ir_genes)

# Exclude genes that are present in PBS_IR_DMXAA and control
exclusive_ure_genes = common_ure_genes.difference(pbs_ir_dmxaa_genes).difference(control_genes)

# Output the results
print("Common genes in Ure_IR and Ure_mock_IR, but not in PBS_IR_DMXAA or control:")
for gene in exclusive_ure_genes:
    print(gene)

# 1. Venn Diagram: Ure conditions vs Non-Ure conditions
ure_conditions = ure_ir_genes.union(ure_mock_ir_genes).union(ure_mock_ir_dmxaa_genes)
non_ure_conditions = pbs_ir_dmxaa_genes.union(control_genes)

plt.figure(figsize=(8, 8))
venn2([ure_conditions, non_ure_conditions],
      set_labels=('Ure Conditions', 'Non-Ure Conditions'),
      set_colors=('skyblue', 'lightgreen'))
plt.title("Venn Diagram: Ure Conditions vs Non-Ure Conditions")
plt.savefig(os.path.join(output_directory, 'venn_ure_vs_non_ure.png'))
plt.show()

# 2. Venn Diagram: Comparing the three Ure conditions
plt.figure(figsize=(10, 10))
venn3([ure_ir_genes, ure_mock_ir_genes, ure_mock_ir_dmxaa_genes],
      set_labels=('Ure_IR', 'Ure_mock_IR', 'Ure_mock_IR_DMXAA'),
      set_colors=('gold', 'lightcoral', 'lightblue'))
plt.title("Venn Diagram: Comparing Ure_IR, Ure_mock_IR, and Ure_mock_IR_DMXAA")
plt.savefig(os.path.join(output_directory, 'venn_ure_conditions_comparison.png'))
plt.show()

# Further filter to keep only the relevant consequences
filtered_df = df[df['SYMBOL'].isin(exclusive_ure_genes)]
filtered_df = filtered_df[filtered_df['Consequence'].isin(['frameshift_variant','missense_variant','frameshift_variant&NMD_transcript_variant'])]

# Generate a color palette with enough distinct colors
unique_consequences = filtered_df['Consequence'].nunique()
palette = sns.color_palette("Set2", unique_consequences)

# Plot 1: Frequency of Consequence Types for Exclusive Ure Genes
plt.figure(figsize=(12, 8))
sns.countplot(y='Consequence', data=filtered_df, order=filtered_df['Consequence'].value_counts().index, palette=palette)
plt.title('Distribution of Consequences for Exclusive Ure Condition Genes')
plt.xlabel('Count')
plt.ylabel('Consequence')
plt.savefig(os.path.join(output_directory, 'exclusive_ure_genes_consequence_distribution.png'))
plt.show()

# Plot 2: Heatmap of Gene Symbols by Consequence
heatmap_data = pd.crosstab(filtered_df['SYMBOL'], filtered_df['Consequence'])
plt.figure(figsize=(14, 12))
sns.heatmap(heatmap_data, cmap="YlGnBu", annot=True, fmt="d")
plt.title('Heatmap of Gene Symbols by Consequence for Exclusive Ure Condition Genes')
plt.xlabel('Consequence')
plt.ylabel('Gene Symbol')
plt.savefig(os.path.join(output_directory, 'exclusive_ure_genes_heatmap.png'))
plt.show()

# Stacked bar plot for frameshift and missense variants
pivot_df = filtered_df.pivot_table(index='SYMBOL', columns='Consequence', aggfunc='size', fill_value=0)

# Stacked bar plot
pivot_df.plot(kind='barh', stacked=True, color=['#66c2a5', '#fc8d62'], figsize=(14, 10))
plt.title('Exclusive Ure Condition Genes with Frameshift/ Missense Variants (Stacked)')
plt.xlabel('Count')
plt.ylabel('Gene Symbol')
plt.legend(title='Consequence', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.savefig(os.path.join(output_directory, 'exclusive_ure_genes_stacked_barplot.png'))
plt.show()

# Heatmap for frameshift and missense variants
plt.figure(figsize=(14, 10))
sns.heatmap(pivot_df, annot=True, cmap="YlGnBu", cbar=True, linewidths=.5)
plt.title('Heatmap of Frameshift/ Missense Variants in Exclusive Ure Condition Genes')
plt.xlabel('Consequence')
plt.ylabel('Gene Symbol')
plt.savefig(os.path.join(output_directory, 'exclusive_ure_genes_stacked_heatmap.png'))
plt.show()

# Dot Plot of Frameshift/ Missense Variants in Exclusive Ure Condition Genes
plt.figure(figsize=(14, 10))
sns.stripplot(x='Consequence', y='SYMBOL', data=filtered_df, size=10, palette="Set2", jitter=True)
plt.title('Dot Plot of Frameshift/ Missense Variants in Exclusive Ure Condition Genes')
plt.xlabel('Consequence')
plt.ylabel('Gene Symbol')
plt.savefig(os.path.join(output_directory, 'exclusive_ure_genes_dotplot.png'))
plt.show()
