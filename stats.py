
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load your CSV file into a pandas DataFrame
df = pd.read_csv('combined_variants_with_sample_info.csv')

# Extract unique genes (CSQ_SYMBOL) for each sample group
ure_ir_genes = set(df[df['Sample_Group'] == 'Ure_IR']['CSQ_SYMBOL'].unique())
ure_mock_ir_genes = set(df[df['Sample_Group'] == 'Ure_mock_IR']['CSQ_SYMBOL'].unique())
ure_mock_ir_dmxaa_genes = set(df[df['Sample_Group'] == 'Ure_mock_IR_DMXAA']['CSQ_SYMBOL'].unique())
pbs_ir_dmxaa_genes = set(df[df['Sample_Group'] == 'PBS_IR_DMXAA']['CSQ_SYMBOL'].unique())
control_genes = set(df[df['Sample_Group'] == 'control']['CSQ_SYMBOL'].unique())

# Find genes that are common in all three Ure conditions
common_ure_genes = ure_ir_genes.intersection(ure_mock_ir_genes).intersection(ure_mock_ir_dmxaa_genes)

# Find which of these common genes are NOT present in PBS_IR_DMXAA and control
exclusive_ure_genes = common_ure_genes.difference(pbs_ir_dmxaa_genes).difference(control_genes)

# Filter the original dataframe to keep only the rows corresponding to the exclusive Ure genes
exclusive_ure_df = df[df['CSQ_SYMBOL'].isin(exclusive_ure_genes)]

# Generate a color palette with enough distinct colors
unique_consequences = exclusive_ure_df['CSQ_Consequence'].nunique()
palette = sns.color_palette("tab20", unique_consequences)

# Plot setup
plt.figure(figsize=(14, 10))

# Count plot using seaborn with adjusted bar width
bars = sns.countplot(y='CSQ_SYMBOL', hue='CSQ_Consequence', data=exclusive_ure_df, 
                     palette=palette, linewidth=0)

# Adjust the width of the bars
for bar in bars.patches:
    bar.set_height(bar.get_height() * 1.5)  # Increase the height by 1.5 times

# Add plot titles and labels
plt.title('Exclusive Ure Condition Genes by Consequence')
plt.xlabel('Count')
plt.ylabel('Gene Symbol')
plt.legend(title='CSQ Consequence', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load your CSV file into a pandas DataFrame
df = pd.read_csv('combined_variants_with_sample_info.csv')

# Extract unique genes (CSQ_SYMBOL) for each sample group
ure_ir_genes = set(df[df['Sample_Group'] == 'Ure_IR']['CSQ_SYMBOL'].unique())
ure_mock_ir_genes = set(df[df['Sample_Group'] == 'Ure_mock_IR']['CSQ_SYMBOL'].unique())
ure_mock_ir_dmxaa_genes = set(df[df['Sample_Group'] == 'Ure_mock_IR_DMXAA']['CSQ_SYMBOL'].unique())
pbs_ir_dmxaa_genes = set(df[df['Sample_Group'] == 'PBS_IR_DMXAA']['CSQ_SYMBOL'].unique())
control_genes = set(df[df['Sample_Group'] == 'control']['CSQ_SYMBOL'].unique())

# Find genes that are common in all three Ure conditions
common_ure_genes = ure_ir_genes.intersection(ure_mock_ir_genes).intersection(ure_mock_ir_dmxaa_genes)

# Find which of these common genes are NOT present in PBS_IR_DMXAA and control
exclusive_ure_genes = common_ure_genes.difference(pbs_ir_dmxaa_genes).difference(control_genes)

# Filter the original dataframe to keep only the rows corresponding to the exclusive Ure genes
exclusive_ure_df = df[df['CSQ_SYMBOL'].isin(exclusive_ure_genes)]

# Further filter to keep only the relevant consequences
filtered_df = exclusive_ure_df[exclusive_ure_df['CSQ_Consequence'].isin(['frameshift_variant', 'frameshift_variant&NMD_transcript_variant'])]

# Generate a color palette with two distinct colors
palette = sns.color_palette("Set2", 2)

# Plot setup
plt.figure(figsize=(14, 10))

# Count plot using seaborn with adjusted bar width
bars = sns.countplot(y='CSQ_SYMBOL', hue='CSQ_Consequence', data=filtered_df, 
                     palette=palette, linewidth=0)

# Adjust the width of the bars
for bar in bars.patches:
    bar.set_height(bar.get_height() * 1.5)  # Increase the height by 1.5 times

# Add plot titles and labels
plt.title('Exclusive Ure Condition Genes with Frameshift Variants')
plt.xlabel('Count')
plt.ylabel('Gene Symbol')
plt.legend(title='CSQ Consequence', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Assuming filtered_df is the DataFrame filtered for frameshift variants
pivot_df = filtered_df.pivot_table(index='CSQ_SYMBOL', columns='CSQ_Consequence', aggfunc='size', fill_value=0)

# Stacked bar plot
pivot_df.plot(kind='barh', stacked=True, color=['#66c2a5', '#fc8d62'], figsize=(14, 10))
plt.title('Exclusive Ure Condition Genes with Frameshift Variants (Stacked)')
plt.xlabel('Count')
plt.ylabel('Gene Symbol')
plt.legend(title='CSQ Consequence', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Assuming filtered_df is the DataFrame filtered for frameshift variants
pivot_df = filtered_df.pivot_table(index='CSQ_SYMBOL', columns='CSQ_Consequence', aggfunc='size', fill_value=0)

# Heatmap
plt.figure(figsize=(14, 10))
sns.heatmap(pivot_df, annot=True, cmap="YlGnBu", cbar=True, linewidths=.5)
plt.title('Heatmap of Frameshift Variants in Exclusive Ure Condition Genes')
plt.xlabel('CSQ Consequence')
plt.ylabel('Gene Symbol')
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Assuming filtered_df is the DataFrame filtered for frameshift variants
plt.figure(figsize=(14, 10))
sns.stripplot(x='CSQ_Consequence', y='CSQ_SYMBOL', data=filtered_df, size=10, palette="Set2", jitter=True)
plt.title('Dot Plot of Frameshift Variants in Exclusive Ure Condition Genes')
plt.xlabel('CSQ Consequence')
plt.ylabel('Gene Symbol')
plt.show()
