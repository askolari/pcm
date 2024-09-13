import pandas as pd
import os

# Define the directory containing the merged files
merged_files_dir = '/mnt/c/Users/agsko/dev/pcm/merged_files/'

# List of all indel files to process
indel_files = [
    "merged_Exome-seq_Sample_1_S3_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_2_S4_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_3_S5_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_4_S6_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_5_S7_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_6_S8_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_7_S9_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_8_S10_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_9_S11_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_10_S12_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_11_S13_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_12_S14_indels_filtered_lifted.csv",
    "merged_Exome-seq_Sample_13_S15_indels_filtered_lifted.csv"
]

# Define sample metadata (mapping Sample ID to Group)
sample_metadata = {
    'S1': 'control',
    'S2': 'control',
    'S3': 'Ure_mock_IR',
    'S4': 'Ure_mock_IR',
    'S5': 'Ure_mock_IR',
    'S6': 'Ure_mock_IR',
    'S7': 'Ure_mock_IR_DMXAA',
    'S8': 'Ure_mock_IR_DMXAA',
    'S9': 'Ure_mock_IR_DMXAA',
    'S10': 'Ure_mock_IR_DMXAA',
    'S11': 'Ure_IR',
    'S12': 'Ure_IR',
    'S13': 'Ure_IR',
    'S14': 'Ure_IR',
    'S15': 'PBS_IR_DMXAA'
}

# Function to extract the first value before a comma
def extract_first_value(val):
    if isinstance(val, str):
        return float(val.split(',')[0])
    return val  # If it's already a number, just return it

# Function to calculate VAF for indels
def calculate_vaf_indel(row):
    ref_counts = row['TUMOR.TAR']
    alt_counts = row['TUMOR.TIR']

    if pd.isna(alt_counts) or (alt_counts + ref_counts) == 0:
        return None
    
    return alt_counts / (alt_counts + ref_counts) * 100

# Initialize an empty DataFrame to hold all the data
combined_indel_df = pd.DataFrame()

# Iterate over each indel file and process it
for file_name in indel_files:
    # Load the specific file into a DataFrame
    file_path = os.path.join(merged_files_dir, file_name)
    df = pd.read_csv(file_path)

    # Apply the extraction to the relevant columns
    df['TUMOR.DP'] = df['TUMOR.DP'].apply(extract_first_value)
    df['TUMOR.TAR'] = df['TUMOR.TAR'].apply(extract_first_value)
    df['TUMOR.TIR'] = df['TUMOR.TIR'].apply(extract_first_value)

    # Apply the function to calculate VAF for each row
    df['VAF'] = df.apply(calculate_vaf_indel, axis=1)

    # Extract Sample_ID from the file name (e.g., "S3" from "merged_Exome-seq_Sample_1_S3_indels_filtered_lifted.csv")
    sample_id = file_name.split('_')[4]

    # Add Sample_ID and Sample_Group to the DataFrame
    df['Sample_ID'] = sample_id
    df['Sample_Group'] = sample_metadata.get(sample_id, 'Unknown')

    # Save the updated DataFrame with the VAF column to a new CSV file
    output_file_name = file_name.replace('.csv', '_with_VAF.csv')
    output_file_path = os.path.join(merged_files_dir, output_file_name)
    df.to_csv(output_file_path, index=False)
    
    print(f"File processed and saved successfully to {output_file_path}")

    # Append the DataFrame to the combined DataFrame
    combined_indel_df = pd.concat([combined_indel_df, df], ignore_index=True)

# Save the combined DataFrame to a new CSV file
output_file_name = 'combined_all_indels_with_VAF_and_Metadata.csv'
output_file_path = os.path.join(merged_files_dir, output_file_name)
combined_indel_df.to_csv(output_file_path, index=False)

print(f"All indel files merged successfully and saved to {output_file_path}")
