import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Specify the directory where you want to save the images
output_directory = '/mnt/c/Users/agsko/dev/pcm/WES_explor/significant_variants'
merged_files_dir = '/mnt/c/Users/agsko/dev/pcm/merged_files/'

# Ensure the output directory exists
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Load the final dataset into a pandas DataFrame
df = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/merged_files/final_variants_mm10.csv')

# Filter for high and moderate impact variants
significant_impact_df = df[df['IMPACT'].isin(['HIGH', 'MODERATE'])]

# Save the filtered high and moderate impact variants to a new CSV file
significant_impact_file_path = os.path.join(merged_files_dir, 'high_moderate_impact_variants.csv')
significant_impact_df.to_csv(significant_impact_file_path, index=False)
print(f"High and moderate impact variants saved to {significant_impact_file_path}")

# Distribution of Variant Types among High and Moderate Impact Variants
sns.countplot(x='Variant_Type', data=significant_impact_df)
plt.title('Distribution of Variant Types among High and Moderate Impact Variants')
plt.savefig(os.path.join(output_directory, 'significant_impact_variant_type_distribution.png'))
plt.show()

# Distribution of Sample Groups among High and Moderate Impact Variants
sns.countplot(x='Sample_Group', data=significant_impact_df)
plt.title('Distribution of Sample Groups among High and Moderate Impact Variants')
plt.savefig(os.path.join(output_directory, 'significant_impact_sample_group_distribution.png'))
plt.show()

# Scatter plot of VAF and variant type among high and moderate impact variants
plt.figure(figsize=(10, 6))
sns.scatterplot(y='Consequence', x='VAF', hue='Sample_Group', data=significant_impact_df)
plt.title('Consequence vs. VAF for High and Moderate Impact Variants')
plt.savefig(os.path.join(output_directory, 'significant_impact_consequence_vs_vaf.png'))
plt.show()

# Box plot for VAF across different Variant Types for high and moderate impact variants
sns.boxplot(x='Variant_Type', y='VAF', data=significant_impact_df)
plt.title('VAF by Variant Type for High and Moderate Impact Variants')
plt.savefig(os.path.join(output_directory, 'significant_impact_vaf_by_variant_type.png'))
plt.show()

# Filter for Ure_IR and Ure_mock_IR comparison
ure_ir_df = significant_impact_df[significant_impact_df['Sample_Group'] == 'Ure_IR']
ure_mock_ir_df = significant_impact_df[significant_impact_df['Sample_Group'] == 'Ure_mock_IR']

# Save these dataframes
ure_ir_file_path = os.path.join(merged_files_dir, 'Ure_IR_significant_variants.csv')
ure_mock_ir_file_path = os.path.join(merged_files_dir, 'Ure_mock_IR_significant_variants.csv')
ure_ir_df.to_csv(ure_ir_file_path, index=False)
ure_mock_ir_df.to_csv(ure_mock_ir_file_path, index=False)
print(f"Ure_IR significant variants saved to {ure_ir_file_path}")
print(f"Ure_mock_IR significant variants saved to {ure_mock_ir_file_path}")

# Compare the distribution of high and moderate impact variants between Ure_IR and Ure_mock_IR
plt.figure(figsize=(12, 6))
sns.countplot(x='IMPACT', hue='Sample_Group', data=significant_impact_df[significant_impact_df['Sample_Group'].isin(['Ure_IR', 'Ure_mock_IR'])], palette="Set2")
plt.title('Significant Variant Impact Distribution between Ure_IR and Ure_mock_IR')
plt.xlabel('Impact')
plt.ylabel('Count')
plt.savefig(os.path.join(output_directory, 'significant_impact_ure_ir_vs_ure_mock_ir.png'))
plt.show()

# Top Genes in Ure_IR vs. Ure_mock_IR for significant variants
top_genes_ure_ir = ure_ir_df['SYMBOL'].value_counts().head(10)
top_genes_ure_mock_ir = ure_mock_ir_df['SYMBOL'].value_counts().head(10)

plt.figure(figsize=(12, 6))
sns.barplot(x=top_genes_ure_ir.values, y=top_genes_ure_ir.index, palette="Blues_d")
plt.title('Top 10 Genes in Ure_IR (Significant Variants)')
plt.xlabel('Count')
plt.ylabel('Gene')
plt.savefig(os.path.join(output_directory, 'top_10_genes_ure_ir_significant.png'))
plt.show()

plt.figure(figsize=(12, 6))
sns.barplot(x=top_genes_ure_mock_ir.values, y=top_genes_ure_mock_ir.index, palette="Reds_d")
plt.title('Top 10 Genes in Ure_mock_IR (Significant Variants)')
plt.xlabel('Count')
plt.ylabel('Gene')
plt.savefig(os.path.join(output_directory, 'top_10_genes_ure_mock_ir_significant.png'))
plt.show()

# Compare specific consequences between Ure_IR and Ure_mock_IR
consequences_comparison_df = pd.crosstab(significant_impact_df[significant_impact_df['Sample_Group'].isin(['Ure_IR', 'Ure_mock_IR'])]['Consequence'], 
                                         significant_impact_df[significant_impact_df['Sample_Group'].isin(['Ure_IR', 'Ure_mock_IR'])]['Sample_Group'])

plt.figure(figsize=(12, 8))
sns.heatmap(consequences_comparison_df, annot=True, cmap="YlGnBu", fmt='d')
plt.title('Consequence Comparison between Ure_IR and Ure_mock_IR (Significant Variants)')
plt.savefig(os.path.join(output_directory, 'consequence_comparison_ure_ir_vs_ure_mock_ir_significant.png'))
plt.show()

# Compare VAF distribution between Ure_IR and Ure_mock_IR for significant variants
plt.figure(figsize=(10, 6))
sns.boxplot(x='Sample_Group', y='VAF', data=significant_impact_df[significant_impact_df['Sample_Group'].isin(['Ure_IR', 'Ure_mock_IR'])], palette="Set3")
plt.title('VAF Distribution between Ure_IR and Ure_mock_IR (Significant Variants)')
plt.savefig(os.path.join(output_directory, 'vaf_distribution_ure_ir_vs_ure_mock_ir_significant.png'))
plt.show()

