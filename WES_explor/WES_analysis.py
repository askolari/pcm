import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Specify the directory where you want to save the images
output_directory = '/mnt/c/Users/agsko/dev/pcm/WES_explor'

# Ensure the output directory exists
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Load your CSV file into a pandas DataFrame
df = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/merged_files/final_variants_mm10.csv')

# Print the column names
print("Column names in the DataFrame:")
for idx, col in enumerate(df.columns):
    print(f"{idx}: {col}")

# Descriptive statistics
print(df.describe())

# Distribution of Variant Types
sns.countplot(x='Variant_Type', data=df)
plt.title('Distribution of Variant Types')
plt.savefig(os.path.join(output_directory, 'variant_type_distribution.png'))
plt.show()

# Distribution of Sample Groups
sns.countplot(x='Sample_Group', data=df)
plt.title('Distribution of Sample Groups')
plt.savefig(os.path.join(output_directory, 'sample_group_distribution.png'))
plt.show()

# Explore specific columns (e.g., Consequence)
sns.countplot(y='Consequence', data=df, order=df['Consequence'].value_counts().index)
plt.title('Distribution of Consequences')
plt.savefig(os.path.join(output_directory, 'consequence_distribution.png'))
plt.show()

# Box plot for VAF across different Variant Types
sns.boxplot(x='Variant_Type', y='VAF', data=df)
plt.title('Variant Allele Frequencies by Variant Type')
plt.savefig(os.path.join(output_directory, 'vaf_by_variant_type.png'))
plt.show()

# Scatter plot of VAF and variant type, colored by Sample Group
plt.figure(figsize=(10, 6))
sns.scatterplot(y='Consequence', x='VAF', hue='Sample_Group', data=df)
plt.title('Consequence vs. VAF, colored by Sample_Group')
plt.savefig(os.path.join(output_directory, 'consequence_vs_vaf.png'))
plt.show()

# Select only numeric columns for correlation matrix
numeric_df = df.select_dtypes(include=['number'])

# Correlation matrix of numeric variables
plt.figure(figsize=(8, 6))
sns.heatmap(numeric_df.corr(), annot=True, fmt=".2f", cmap='coolwarm')
plt.title('Correlation Matrix')
plt.savefig(os.path.join(output_directory, 'correlation_matrix.png'))
plt.show()

# Focus on specific columns
df = df[['Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Sample_ID', 'Sample_Group', 'Variant_Type']]

# 1. Distribution of Consequence with color differentiation
plt.figure(figsize=(12, 8))
sns.countplot(y='Consequence', data=df, order=df['Consequence'].value_counts().index, palette="coolwarm")
plt.title('Distribution of Consequences')
plt.xlabel('Count')
plt.ylabel('Consequence')
plt.savefig(os.path.join(output_directory, 'consequence_distribution_detailed.png'))
plt.show()

# 2. Distribution of CSQ Impact with color differentiation
plt.figure(figsize=(10, 6))
sns.countplot(x='IMPACT', data=df, palette="Set2")
plt.title('Distribution of Impact')
plt.xlabel('Impact')
plt.ylabel('Count')
plt.savefig(os.path.join(output_directory, 'impact_distribution.png'))
plt.show()

# 3. Count Plot of Variant Types by Sample Group
plt.figure(figsize=(12, 6))
sns.countplot(x='Sample_Group', hue='Variant_Type', data=df, palette="husl")
plt.title('Variant Types by Sample Group')
plt.xlabel('Sample Group')
plt.ylabel('Count')
plt.xticks(rotation=45)
plt.savefig(os.path.join(output_directory, 'variant_types_by_sample_group.png'))
plt.show()

# 4. Top 10 Most Frequent Genes in the Data
top_genes = df['Gene'].value_counts().head(10)
plt.figure(figsize=(12, 6))
sns.barplot(x=top_genes.values, y=top_genes.index, palette="magma")
plt.title('Top 10 Most Frequent Genes')
plt.xlabel('Count')
plt.ylabel('Gene')
plt.savefig(os.path.join(output_directory, 'top_10_genes.png'))
plt.show()

# 5. Relationship between CSQ Impact and Variant Type
plt.figure(figsize=(12, 6))
sns.countplot(x='IMPACT', hue='Variant_Type', data=df, palette="coolwarm")
plt.title('Impact by Variant Type')
plt.xlabel('Impact')
plt.ylabel('Count')
plt.savefig(os.path.join(output_directory, 'impact_by_variant_type.png'))
plt.show()

# 6. Relationship between CSQ Consequence and CSQ Impact
plt.figure(figsize=(14, 10))
sns.heatmap(pd.crosstab(df['Consequence'], df['IMPACT']), cmap="YlGnBu", annot=True, fmt='d')
plt.title('Heatmap of Consequence vs. Impact')
plt.xlabel('Impact')
plt.ylabel('Consequence')
plt.savefig(os.path.join(output_directory, 'consequence_vs_impact_heatmap.png'))
plt.show()

# 7. Top 10 Most Frequent Symbols in the Data
top_symbols = df['SYMBOL'].value_counts().head(10)
plt.figure(figsize=(12, 6))
sns.barplot(x=top_symbols.values, y=top_symbols.index, palette="magma")
plt.title('Top 10 Most Frequent SYMBOLs')
plt.xlabel('Count')
plt.ylabel('SYMBOL')
plt.savefig(os.path.join(output_directory, 'top_10_symbols.png'))
plt.show()


from matplotlib_venn import venn2, venn3

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
