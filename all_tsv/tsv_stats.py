import pandas as pd

# Load the two datasets
df1 = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/all_tsv/combined_tsv_variants.csv')
df2 = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/all_tsv/combined_variants_with_sample_info.csv')

# Identify the common columns for merging
common_columns = ['CHROM', 'POS', 'REF', 'ALT']

# Merge the datasets on the common columns
merged_df = pd.merge(df1, df2, on=common_columns, how='inner')

# Save the merged dataset to a new CSV file
merged_df.to_csv('/mnt/c/Users/agsko/dev/pcm/all_tsv/merged_variants_dataset.csv', index=False)

print("Datasets merged successfully and saved as 'merged_variants_dataset.csv'.")
