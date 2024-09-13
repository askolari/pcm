import pandas as pd
import matplotlib.pyplot as plt
import seaborn as snsfilter

# Load the final dataset
final_dataset_path = '/mnt/c/Users/agsko/dev/pcm/merged_files/final_variants_mm10.csv'  # Replace with your actual file path
df = pd.read_csv(final_dataset_path)

# Filter for frameshift variants
frameshift_df = df[df['Consequence'].str.contains('frameshift_variant', case=False, na=False)]

# Plot 1: Dot plot of frameshift genes by group with VAF on x-axis
plt.figure(figsize=(12, 8))
sns.stripplot(x='VAF', y='SYMBOL', hue='Sample_Group', data=frameshift_df, jitter=True, dodge=True, marker='o', palette="Set2")
plt.title('Frameshift Genes by Group with %VAF')
plt.xlabel('%VAF')
plt.ylabel('Gene')
plt.legend(title='Group', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Save the plot
plot_path = '/mnt/c/Users/agsko/dev/pcm/WES_explor/frameshift_genes_by_group_vaf_dotplot.png'  # Replace with your desired output path
plt.savefig(plot_path)
plt.show()

print(f"Dot plot of frameshift genes by group saved to {plot_path}")

# Plot 1: Bar plot of frameshift genes by group with VAF on x-axis
plt.figure(figsize=(12, 8))
sns.barplot(x='VAF', y='SYMBOL', hue='Sample_Group', data=frameshift_df, palette="Set2", ci=None)
plt.title('Frameshift Genes by Group with %VAF')
plt.xlabel('%VAF')
plt.ylabel('Gene')
plt.legend(title='Group', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Save the plot
plot_path = '/mnt/c/Users/agsko/dev/pcm/WES_explor/frameshift_genes_by_group_vaf_barplot.png'  # Replace with your desired output path
plt.savefig(plot_path)
plt.show()


