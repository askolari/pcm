import pandas as pd
import os

# Define the directory containing the merged files
merged_files_dir = '/mnt/c/Users/agsko/dev/pcm/merged_files/'

# List of all SNV files to process
snv_files = [
    "merged_Exome-seq_Control_1_S1_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Control_2_S2_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_1_S3_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_2_S4_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_3_S5_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_4_S6_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_5_S7_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_6_S8_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_7_S9_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_8_S10_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_9_S11_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_10_S12_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_11_S13_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_12_S14_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Sample_13_S15_snvs_filtered_lifted.csv"
]

# Function to extract the first value before a comma
def extract_first_value(val):
    if isinstance(val, str):
        return float(val.split(',')[0])
    return val  # If it's already a number, just return it


# Function to calculate VAF for SNVs
def calculate_vaf_snv(row):
    if row['ALT_x'] == 'A':
        alt_counts = row['TUMOR.AU']
        ref_counts = row['TUMOR.TU'] + row['TUMOR.GU'] + row['TUMOR.CU']
    elif row['ALT_x'] == 'C':
        alt_counts = row['TUMOR.CU']
        ref_counts = row['TUMOR.TU'] + row['TUMOR.GU'] + row['TUMOR.AU']
    elif row['ALT_x'] == 'G':
        alt_counts = row['TUMOR.GU']
        ref_counts = row['TUMOR.TU'] + row['TUMOR.AU'] + row['TUMOR.CU']
    elif row['ALT_x'] == 'T':
        alt_counts = row['TUMOR.TU']
        ref_counts = row['TUMOR.GU'] + row['TUMOR.AU'] + row['TUMOR.CU']
    else:
        return None

    if pd.isna(alt_counts) or (alt_counts + ref_counts) == 0:
        return None
    
    return alt_counts / (alt_counts + ref_counts) * 100
    

# Iterate over each SNV file
for file_name in snv_files:
    # Load the specific file into a DataFrame
    file_path = os.path.join(merged_files_dir, file_name)
    df = pd.read_csv(file_path)

    # Apply the extraction to the relevant columns
    numeric_columns = ['TUMOR.DP', 'TUMOR.AU', 'TUMOR.CU', 'TUMOR.GU', 'TUMOR.TU']
    df[numeric_columns] = df[numeric_columns].applymap(extract_first_value)

    # Apply the function to calculate VAF for each row
    df['VAF'] = df.apply(calculate_vaf_snv, axis=1)

    # Save the updated DataFrame with the VAF column to a new CSV file
    output_file_name = file_name.replace('.csv', '_with_VAF.csv')
    output_file_path = os.path.join(merged_files_dir, output_file_name)
    df.to_csv(output_file_path, index=False)

    print(f"File processed and saved successfully to {output_file_path}")

import pandas as pd
import os

# Define the directory containing the merged files with VAF
merged_files_dir = '/mnt/c/Users/agsko/dev/pcm/merged_files/'

# List of all SNV files with VAF to merge
snv_files_with_vaf = [
    "merged_Exome-seq_Control_2_S2_snvs_filtered_lifted_with_VAF.csv",
    "merged_Exome-seq_Sample_1_S3_snvs_filtered_lifted_with_VAF.csv",
    "merged_Exome-seq_Sample_2_S4_snvs_filtered_lifted_with_VAF.csv",
    "merged_Exome-seq_Sample_3_S5_snvs_filtered_lifted_with_VAF.csv",
    "merged_Exome-seq_Sample_4_S6_snvs_filtered_lifted_with_VAF.csv",
    "merged_Exome-seq_Sample_5_S7_snvs_filtered_lifted_with_VAF.csv",
    "merged_Exome-seq_Sample_6_S8_snvs_filtered_lifted_with_VAF.csv",
    "merged_Exome-seq_Sample_7_S9_snvs_filtered_lifted_with_VAF.csv",
    "merged_Exome-seq_Sample_8_S10_snvs_filtered_lifted_with_VAF.csv",
    "merged_Exome-seq_Sample_9_S11_snvs_filtered_lifted_with_VAF.csv",
    "merged_Exome-seq_Sample_10_S12_snvs_filtered_lifted_with_VAF.csv",
    "merged_Exome-seq_Sample_11_S13_snvs_filtered_lifted_with_VAF.csv",
    "merged_Exome-seq_Sample_12_S14_snvs_filtered_lifted_with_VAF.csv",
    "merged_Exome-seq_Sample_13_S15_snvs_filtered_lifted_with_VAF.csv"
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

# Initialize an empty DataFrame to hold all the data
combined_snv_df = pd.DataFrame()

# Iterate over each SNV file and add Sample_ID and Sample_Group
for file_name in snv_files_with_vaf:
    file_path = os.path.join(merged_files_dir, file_name)
    
    # Read the current file into a DataFrame
    df = pd.read_csv(file_path)
    
    # Extract Sample_ID from the file name (e.g., "S3" from "merged_Exome-seq_Sample_1_S3_snvs_filtered_lifted_with_VAF.csv")
    sample_id = file_name.split('_')[4]
    
    # Add Sample_ID and Sample_Group to the DataFrame
    df['Sample_ID'] = sample_id
    df['Sample_Group'] = sample_metadata.get(sample_id, 'Unknown')
    
    # Append the DataFrame to the combined DataFrame
    combined_snv_df = pd.concat([combined_snv_df, df], ignore_index=True)

# Save the combined DataFrame to a new CSV file
output_file_name = 'combined_all_snv_with_VAF_and_Metadata.csv'
output_file_path = os.path.join(merged_files_dir, output_file_name)
combined_snv_df.to_csv(output_file_path, index=False)

print(f"All SNV files merged successfully and saved to {output_file_path}")

