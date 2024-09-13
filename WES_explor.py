import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Specify the directory where you want to save the images
output_directory = '/mnt/c/Users/agsko/dev/pcm/WES_explor'

# Specify the file name and path for saving the plot/image
file_path = os.path.join(output_directory, "my_image.png")

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
plt.show()

# Distribution of Sample Groups
sns.countplot(x='Sample_Group', data=df)
plt.title('Distribution of Sample Groups')
plt.show()

# Explore specific columns (e.g., CSQ_Consequence)
sns.countplot(y='Consequence', data=df, order=df['Consequence'].value_counts().index)
plt.title('Distribution of Consequences')
plt.show()

# Box plot for QUAL scores across different Variant Types
sns.boxplot(x='Variant_Type', y='QUAL', data=df)
plt.title('QUAL Scores by Variant Type')
plt.show()

# Scatter plot of Position (POS) vs. QUAL
plt.figure(figsize=(10, 6))
sns.scatterplot(x='POS', y='QUAL', hue='Variant_Type', data=df)
plt.title('Position vs. QUAL, colored by Variant Type')
plt.show()

# Correlation matrix of numeric variables
plt.figure(figsize=(8, 6))
sns.heatmap(df.corr(), annot=True, fmt=".2f")
plt.title('Correlation Matrix')
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load your CSV file into a pandas DataFrame
df = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/merged_files/final_variants_mm10.csv')

# Focus on the specified columns
df = df[['Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Sample_ID', 'Sample_Group', 'Variant_Type']]

# Set the Seaborn style for the plots
sns.set(style="whitegrid")

# 1. Distribution of Consequence with color differentiation
plt.figure(figsize=(12, 8))
sns.countplot(y='Consequence', data=df, order=df['Consequence'].value_counts().index, palette="coolwarm")
plt.title('Distribution of Consequences')
plt.xlabel('Count')
plt.ylabel('Consequence')
plt.show()

# 2. Distribution of CSQ Impact with color differentiation
plt.figure(figsize=(10, 6))
sns.countplot(x='IMPACT', data=df, palette="Set2")
plt.title('Distribution of Impact')
plt.xlabel('Impact')
plt.ylabel('Count')
plt.show()

# 3. Count Plot of Variant Types by Sample Group
plt.figure(figsize=(12, 6))
sns.countplot(x='Sample_Group', hue='Variant_Type', data=df, palette="husl")
plt.title('Variant Types by Sample Group')
plt.xlabel('Sample Group')
plt.ylabel('Count')
plt.xticks(rotation=45)
plt.show()

# 4. Top 10 Most Frequent Genes in the Data
top_genes = df['Gene'].value_counts().head(10)
plt.figure(figsize=(12, 6))
sns.barplot(x=top_genes.values, y=top_genes.index, palette="magma")
plt.title('Top 10 Most Frequent Genes')
plt.xlabel('Count')
plt.ylabel('Gene')
plt.show()

# 5. Relationship between CSQ Impact and Variant Type
plt.figure(figsize=(12, 6))
sns.countplot(x='IMPACT', hue='Variant_Type', data=df, palette="coolwarm")
plt.title('Impact by Variant Type')
plt.xlabel('Impact')
plt.ylabel('Count')
plt.show()

# 6. Relationship between CSQ Consequence and CSQ Impact
plt.figure(figsize=(14, 10))
sns.heatmap(pd.crosstab(df['Consequence'], df['IMPACT']), cmap="YlGnBu", annot=True, fmt='d')
plt.title('Heatmap of CSQ Consequence vs. CSQ Impact')
plt.xlabel('Impact')
plt.ylabel('Consequence')
plt.show()

# Top 10 Most Frequent CSQ_SYMBOLs in the Data
top_symbols = df['SYMBOL'].value_counts().head(10)
plt.figure(figsize=(12, 6))
sns.barplot(x=top_symbols.values, y=top_symbols.index, palette="magma")
plt.title('Top 10 Most Frequent SYMBOLs')
plt.xlabel('Count')
plt.ylabel('SYMBOL')
plt.show()


pip install matplotlib-venn

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load your CSV file into a pandas DataFrame
df = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/merged_files/final_variants_mm10.csv')

# Focus on the specified columns
df = df[['Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Sample_ID', 'Sample_Group', 'Variant_Type', 'VAF']]

# Set the Seaborn style for the plots
sns.set(style="whitegrid")

import pandas as pd

# Load your CSV file into a pandas DataFrame
df = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/merged_files/final_variants_mm10.csv')

# Focus on the relevant columns
df = df[['SYMBOL', 'Sample_Group']]

# Split the DataFrame by Sample Group
ure_ir_genes = df[df['Sample_Group'] == 'Ure_IR']['SYMBOL'].unique()
ure_mock_ir_genes = df[df['Sample_Group'] == 'Ure_mock_IR']['SYMBOL'].unique()
pbs_ir_dmxaa_genes = df[df['Sample_Group'] == 'PBS_IR_DMXAA']['SYMBOL'].unique()
control_genes = df[df['Sample_Group'] == 'control']['SYMBOL'].unique()

# Find common genes in Ure_IR and Ure_mock_IR
common_ure_genes = set(ure_ir_genes).intersection(set(ure_mock_ir_genes))

# Exclude genes that are present in PBS_IR_DMXAA and control
exclusive_ure_genes = common_ure_genes.difference(set(pbs_ir_dmxaa_genes)).difference(set(control_genes))

# Output the results
print("Common genes in Ure_IR and Ure_mock_IR, but not in PBS_IR_DMXAA or control:")
for gene in exclusive_ure_genes:
    print(gene)

from matplotlib_venn import venn2, venn3
import matplotlib.pyplot as plt

# Load your CSV file into a pandas DataFrame
df = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/merged_files/final_variants_mm10.csv')

# Extract unique genes (CSQ_SYMBOL) for each sample group
ure_ir_genes = set(df[df['Sample_Group'] == 'Ure_IR']['SYMBOL'].unique())
ure_mock_ir_genes = set(df[df['Sample_Group'] == 'Ure_mock_IR']['SYMBOL'].unique())
ure_mock_ir_dmxaa_genes = set(df[df['Sample_Group'] == 'Ure_mock_IR_DMXAA']['SYMBOL'].unique())
pbs_ir_dmxaa_genes = set(df[df['Sample_Group'] == 'PBS_IR_DMXAA']['SYMBOL'].unique())
control_genes = set(df[df['Sample_Group'] == 'control']['SYMBOL'].unique())

# 1. Venn Diagram: Ure conditions vs Non-Ure conditions
ure_conditions = ure_ir_genes.union(ure_mock_ir_genes).union(ure_mock_ir_dmxaa_genes)
non_ure_conditions = pbs_ir_dmxaa_genes.union(control_genes)

plt.figure(figsize=(8, 8))
venn2([ure_conditions, non_ure_conditions],
      set_labels=('Ure Conditions', 'Non-Ure Conditions'),
      set_colors=('skyblue', 'lightgreen'))
plt.title("Venn Diagram: Ure Conditions vs Non-Ure Conditions")
plt.show()

# 2. Venn Diagram: Comparing the three Ure conditions
plt.figure(figsize=(10, 10))
venn3([ure_ir_genes, ure_mock_ir_genes, ure_mock_ir_dmxaa_genes],
      set_labels=('Ure_IR', 'Ure_mock_IR', 'Ure_mock_IR_DMXAA'),
      set_colors=('gold', 'lightcoral', 'lightblue'))
plt.title("Venn Diagram: Comparing Ure_IR, Ure_mock_IR, and Ure_mock_IR_DMXAA")
plt.show()


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load your CSV file into a pandas DataFrame
df = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/merged_files/final_variants_mm10.csv')

# Focus on the specified columns
df = df[['Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Sample_ID', 'Sample_Group', 'Variant_Type', 'VAF']]

# Set the Seaborn style for the plots
sns.set(style="whitegrid")

import pandas as pd

# Load your CSV file into a pandas DataFrame
df = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/merged_files/final_variants_mm10.csv')

# Focus on the relevant columns
df = df[['SYMBOL', 'Sample_Group']]

# Split the DataFrame by Sample Group
ure_ir_genes = df[df['Sample_Group'] == 'Ure_IR']['SYMBOL'].unique()
ure_mock_ir_genes = df[df['Sample_Group'] == 'Ure_mock_IR']['SYMBOL'].unique()
pbs_ir_dmxaa_genes = df[df['Sample_Group'] == 'PBS_IR_DMXAA']['SYMBOL'].unique()
control_genes = df[df['Sample_Group'] == 'control']['SYMBOL'].unique()

import pandas as pd

# Load your CSV file into a pandas DataFrame
df = pd.read_csv('combined_variants_with_sample_info.csv')

# Extract unique genes (CSQ_SYMBOL) for each sample group
ure_ir_genes = set(df[df['Sample_Group'] == 'Ure_IR']['SYMBOL'].unique())
ure_mock_ir_genes = set(df[df['Sample_Group'] == 'Ure_mock_IR']['SYMBOL'].unique())
ure_mock_ir_dmxaa_genes = set(df[df['Sample_Group'] == 'Ure_mock_IR_DMXAA']['SYMBOL'].unique())
pbs_ir_dmxaa_genes = set(df[df['Sample_Group'] == 'PBS_IR_DMXAA']['SYMBOL'].unique())
control_genes = set(df[df['Sample_Group'] == 'control']['SYMBOL'].unique())

# Find genes that are common in all three Ure conditions
common_ure_genes = ure_ir_genes.intersection(ure_mock_ir_genes).intersection(ure_mock_ir_dmxaa_genes)

# Find which of these common genes are NOT present in PBS_IR_DMXAA and control
exclusive_ure_genes = common_ure_genes.difference(pbs_ir_dmxaa_genes).difference(control_genes)

# Output the results
print(f"Number of genes common in all Ure conditions: {len(common_ure_genes)}")
print(f"Number of these genes NOT present in PBS and control: {len(exclusive_ure_genes)}")
print("Exclusive Ure condition genes:")
for gene in exclusive_ure_genes:
    print(gene)

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load your CSV file into a pandas DataFrame
df = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/merged_files/final_variants_mm10.csv')

# Extract unique genes (CSQ_SYMBOL) for each sample group
ure_ir_genes = set(df[df['Sample_Group'] == 'Ure_IR']['SYMBOL'].unique())
ure_mock_ir_genes = set(df[df['Sample_Group'] == 'Ure_mock_IR']['SYMBOL'].unique())
ure_mock_ir_dmxaa_genes = set(df[df['Sample_Group'] == 'Ure_mock_IR_DMXAA']['SYMBOL'].unique())
pbs_ir_dmxaa_genes = set(df[df['Sample_Group'] == 'PBS_IR_DMXAA']['SYMBOL'].unique())
control_genes = set(df[df['Sample_Group'] == 'control']['SYMBOL'].unique())

# Find genes that are common in all three Ure conditions
common_ure_genes = ure_ir_genes.intersection(ure_mock_ir_genes).intersection(ure_mock_ir_dmxaa_genes)

# Find which of these common genes are NOT present in PBS_IR_DMXAA and control
exclusive_ure_genes = common_ure_genes.difference(pbs_ir_dmxaa_genes).difference(control_genes)


# Filter the original dataframe to keep only the rows corresponding to the exclusive Ure genes
exclusive_ure_df = df[df['SYMBOL'].isin(exclusive_ure_genes)]

# Plot 1: Frequency of Consequence Types for Exclusive Ure Genes
plt.figure(figsize=(12, 8))
sns.countplot(y='Consequence', data=exclusive_ure_df, order=exclusive_ure_df['Consequence'].value_counts().index, palette="coolwarm")
plt.title('Distribution of Consequences for Exclusive Ure Condition Genes')
plt.xlabel('Count')
plt.ylabel('Consequence')
plt.show()

# Plot 3: Heatmap of Gene Symbols by Consequence
heatmap_data = pd.crosstab(exclusive_ure_df['SYMBOL'], exclusive_ure_df['Consequence'])
plt.figure(figsize=(14, 12))
sns.heatmap(heatmap_data, cmap="YlGnBu", annot=True, fmt="d")
plt.title('Heatmap of Gene Symbols by Consequence for Exclusive Ure Condition Genes')
plt.xlabel('Consequence')
plt.ylabel('Gene Symbol')
plt.show()


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load your CSV file into a pandas DataFrame
df = pd.read_csv('combined_variants_with_sample_info.csv')

# Extract unique genes (CSQ_SYMBOL) for each sample group
ure_ir_genes = set(df[df['Sample_Group'] == 'Ure_IR']['SYMBOL'].unique())
ure_mock_ir_genes = set(df[df['Sample_Group'] == 'Ure_mock_IR']['SYMBOL'].unique())
ure_mock_ir_dmxaa_genes = set(df[df['Sample_Group'] == 'Ure_mock_IR_DMXAA']['SYMBOL'].unique())
pbs_ir_dmxaa_genes = set(df[df['Sample_Group'] == 'PBS_IR_DMXAA']['SYMBOL'].unique())
control_genes = set(df[df['Sample_Group'] == 'control']['SYMBOL'].unique())

# Find genes that are common in all three Ure conditions
common_ure_genes = ure_ir_genes.intersection(ure_mock_ir_genes).intersection(ure_mock_ir_dmxaa_genes)

# Find which of these common genes are NOT present in PBS_IR_DMXAA and control
exclusive_ure_genes = common_ure_genes.difference(pbs_ir_dmxaa_genes).difference(control_genes)

# Filter the original dataframe to keep only the rows corresponding to the exclusive Ure genes
exclusive_ure_df = df[df['SYMBOL'].isin(exclusive_ure_genes)]

# Generate a color palette with enough distinct colors
unique_consequences = exclusive_ure_df['Consequence'].nunique()
palette = sns.color_palette("tab20", unique_consequences)

# Plot setup
plt.figure(figsize=(14, 10))

# Count plot using seaborn with adjusted bar width
bars = sns.countplot(y='SYMBOL', hue='Consequence', data=exclusive_ure_df, 
                     palette=palette, linewidth=0)

# Adjust the width of the bars
for bar in bars.patches:
    bar.set_height(bar.get_height() * 1.5)  # Increase the height by 1.5 times

# Add plot titles and labels
plt.title('Exclusive Ure Condition Genes by Consequence')
plt.xlabel('Count')
plt.ylabel('Gene Symbol')
plt.legend(title='Consequence', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load your CSV file into a pandas DataFrame
df = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/merged_files/final_variants_mm10.csv')

# Extract unique genes (CSQ_SYMBOL) for each sample group
ure_ir_genes = set(df[df['Sample_Group'] == 'Ure_IR']['SYMBOL'].unique())
ure_mock_ir_genes = set(df[df['Sample_Group'] == 'Ure_mock_IR']['SYMBOL'].unique())
ure_mock_ir_dmxaa_genes = set(df[df['Sample_Group'] == 'Ure_mock_IR_DMXAA']['SYMBOL'].unique())
pbs_ir_dmxaa_genes = set(df[df['Sample_Group'] == 'PBS_IR_DMXAA']['SYMBOL'].unique())
control_genes = set(df[df['Sample_Group'] == 'control']['SYMBOL'].unique())

# Find genes that are common in all three Ure conditions
common_ure_genes = ure_ir_genes.intersection(ure_mock_ir_genes).intersection(ure_mock_ir_dmxaa_genes)

# Find which of these common genes are NOT present in PBS_IR_DMXAA and control
exclusive_ure_genes = common_ure_genes.difference(pbs_ir_dmxaa_genes).difference(control_genes)

# Filter the original dataframe to keep only the rows corresponding to the exclusive Ure genes
exclusive_ure_df = df[df['SYMBOL'].isin(exclusive_ure_genes)]

# Further filter to keep only the relevant consequences
filtered_df = exclusive_ure_df[exclusive_ure_df['Consequence'].isin(['frameshift_variant','missense_variant','frameshift_variant&NMD_transcript_variant'])]

# Generate a color palette with two distinct colors
palette = sns.color_palette("Set2", 2)

# Plot setup
plt.figure(figsize=(14, 10))

# Count plot using seaborn with adjusted bar width
bars = sns.countplot(y='SYMBOL', hue='Consequence', data=filtered_df, 
                     palette=palette, linewidth=0)

# Adjust the width of the bars
for bar in bars.patches:
    bar.set_height(bar.get_height() * 1.5)  # Increase the height by 1.5 times

# Add plot titles and labels
plt.title('Exclusive Ure Condition Genes with Frameshift/ Missense Variants')
plt.xlabel('Count')
plt.ylabel('Gene Symbol')
plt.legend(title='Consequence', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Assuming filtered_df is the DataFrame filtered for frameshift variants
pivot_df = filtered_df.pivot_table(index='SYMBOL', columns='Consequence', aggfunc='size', fill_value=0)

# Stacked bar plot
pivot_df.plot(kind='barh', stacked=True, color=['#66c2a5', '#fc8d62'], figsize=(14, 10))
plt.title('Exclusive Ure Condition Genes with Frameshift/ Missense Variants (Stacked)')
plt.xlabel('Count')
plt.ylabel('Gene Symbol')
plt.legend(title='Consequence', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Assuming filtered_df is the DataFrame filtered for frameshift variants
pivot_df = filtered_df.pivot_table(index='SYMBOL', columns='Consequence', aggfunc='size', fill_value=0)

# Heatmap
plt.figure(figsize=(14, 10))
sns.heatmap(pivot_df, annot=True, cmap="YlGnBu", cbar=True, linewidths=.5)
plt.title('Heatmap of Frameshift/ Missense Variants in Exclusive Ure Condition Genes')
plt.xlabel('Consequence')
plt.ylabel('Gene Symbol')
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Assuming filtered_df is the DataFrame filtered for frameshift variants
plt.figure(figsize=(14, 10))
sns.stripplot(x='Consequence', y='SYMBOL', data=filtered_df, size=10, palette="Set2", jitter=True)
plt.title('Dot Plot of Frameshift/ Missense Variants in Exclusive Ure Condition Genes')
plt.xlabel('Consequence')
plt.ylabel('Gene Symbol')
plt.show()

