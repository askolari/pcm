import pandas as pd
import os

# Define the directory containing the merged files
merged_files_dir = '/mnt/c/Users/agsko/dev/pcm/merged_files/'

# List of all merged CSV files
merged_files = [
    "merged_Exome-seq_Control_2_S2_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_1_S3_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_1_S3_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_2_S4_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_2_S4_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_3_S5_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_3_S5_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_4_S6_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_4_S6_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_5_S7_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_5_S7_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_6_S8_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_6_S8_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_7_S9_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_7_S9_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_8_S10_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_8_S10_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_9_S11_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_9_S11_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_10_S12_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_10_S12_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_11_S13_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_11_S13_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_12_S14_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_12_S14_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_13_S15_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_13_S15_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Control_1_S1_snvs_filtered_lifted.csv"
]

# Initialize an empty DataFrame to hold all data
combined_df = pd.DataFrame()

# Iterate through each file and process it
for file_name in merged_files:
    file_path = os.path.join(merged_files_dir, file_name)
    
    # Read the CSV file
    df = pd.read_csv(file_path)
    
    # Extract Sample_ID and Variant_Type from the file name
    parts = file_name.split('_')
    sample_id = parts[3]  # Sample_ID is the 4th element, e.g., S1, S2, etc.
    variant_type = 'snv' if 'snvs' in file_name else 'indel'  # Determine variant type based on file name
    
    # Add these as new columns to the DataFrame
    df['Sample_ID'] = sample_id
    df['Variant_Type'] = variant_type
    
    # Append to the combined DataFrame
    combined_df = pd.concat([combined_df, df], ignore_index=True)

# Save the combined DataFrame to a new CSV file
output_file_path = '/mnt/c/Users/agsko/dev/pcm/combined_all_samples.csv'
combined_df.to_csv(output_file_path, index=False)

print(f"Combined dataset saved as: {output_file_path}")
