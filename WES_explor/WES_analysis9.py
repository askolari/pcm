import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from statsmodels.stats.contingency_tables import Table2x2
import numpy as np

# Load your dataset
df = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/merged_files/final_variants_mm10.csv')

# Filter to include genes with the most mutations
top_genes = df['SYMBOL'].value_counts().head(20).index
filtered_df = df[df['SYMBOL'].isin(top_genes)]

# 1. Bar Plot: Variant Count by Gene
plt.figure(figsize=(12, 8))
sns.countplot(y='SYMBOL', data=filtered_df, palette="coolwarm", order=filtered_df['SYMBOL'].value_counts().index)
plt.title('Top 20 Genes by Mutation Count')
plt.xlabel('Mutation Count')
plt.ylabel('Gene')
plt.xticks(rotation=45)
plt.show()

# 2. Simplified Heatmap
heatmap_data = pd.crosstab(filtered_df['SYMBOL'], filtered_df['Sample_Group'])
plt.figure(figsize=(12, 8))
sns.heatmap(heatmap_data, annot=True, cmap="YlGnBu", cbar=True, linewidths=.5)
plt.title('Heatmap of Gene Mutations Across Sample Groups')
plt.xlabel('Sample Group')
plt.ylabel('Gene Symbol')
plt.show()


