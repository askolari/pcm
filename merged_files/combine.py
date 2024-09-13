import pandas as pd
import os

# Define the directory containing the merged files
merged_files_dir = '/mnt/c/Users/agsko/dev/pcm/merged_files/'

# Paths to the combined SNV and indel files
snv_file_path = os.path.join(merged_files_dir, 'combined_all_snv_with_VAF_and_Metadata.csv')
indel_file_path = os.path.join(merged_files_dir, 'combined_all_indels_with_VAF_and_Metadata.csv')

# Load the SNV and indel datasets
snv_df = pd.read_csv(snv_file_path)
indel_df = pd.read_csv(indel_file_path)

# Add the Variant_Type column
snv_df['Variant_Type'] = 'SNV'
indel_df['Variant_Type'] = 'INDEL'

# Combine the two datasets
final_combined_df = pd.concat([snv_df, indel_df], ignore_index=True)

# Save the combined dataset to a new CSV file
output_file_name = 'final_variants_mm10.csv'
output_file_path = os.path.join(merged_files_dir, output_file_name)
final_combined_df.to_csv(output_file_path, index=False)

print(f"Final combined variants dataset saved successfully to {output_file_path}")
