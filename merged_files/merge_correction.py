import pandas as pd
import os

# Define the directory containing the merged files
merged_files_dir = '/mnt/c/Users/agsko/dev/pcm/merged_files/'

# List of indel and snv files
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

snv_files = [
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
    "merged_Exome-seq_Sample_13_S15_snvs_filtered_lifted.csv",
    "merged_Exome-seq_Control_1_S1_snvs_filtered_lifted.csv"
]

# Columns to drop (adjust as needed)
columns_to_drop = ['MANE_SELECT', 'MANE_PLUS_CLINICAL']

# Function to calculate VAF for indels
def calculate_vaf_indel(row):
    alt_reads = row['TUMOR.TAR']  # Assuming TUMOR.TAR is the read supporting the indel
    if pd.isna(alt_reads) or row['TUMOR.DP'] == 0 or pd.isna(row['TUMOR.DP']):
        return None
    return alt_reads / row['TUMOR.DP']

def process_indel_files(file_list, output_file_name):
    combined_df = pd.DataFrame()

    for file_name in file_list:
        file_path = os.path.join(merged_files_dir, file_name)
        
        # Read the CSV file
        df = pd.read_csv(file_path)
        
        # Drop unnecessary columns
        df = df.drop(columns=columns_to_drop, errors='ignore')
        
        # Convert relevant columns to numeric
        numeric_columns = ['TUMOR.DP', 'TUMOR.TAR']
        df[numeric_columns] = df[numeric_columns].apply(pd.to_numeric, errors='coerce')
        
        # Calculate VAF for indels
        df['VAF'] = df.apply(calculate_vaf_indel, axis=1)
        
        # Append to the combined DataFrame
        combined_df = pd.concat([combined_df, df], ignore_index=True)

    # Save the combined DataFrame to a new CSV file
    output_file_path = os.path.join(merged_files_dir, output_file_name)
    combined_df.to_csv(output_file_path, index=False)
    print(f"Combined indels dataset saved as: {output_file_path}")

# Process all indel files
process_indel_files(indel_files, 'all_indels_dataset.csv')