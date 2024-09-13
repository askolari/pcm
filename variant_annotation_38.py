import pandas as pd
import os

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

# Function to add sample metadata
def add_sample_metadata(file_path):
    # Extract the sample ID from the file name
    sample_id = os.path.basename(file_path).split('_')[2]
    group = sample_metadata.get(sample_id, 'unknown')
    
    return sample_id, group

# List of parsed CSV file paths
parsed_files = [
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_12_S14_indels_filtered_lifted_annotated.vep_parsed.csv",
    # Add the other parsed CSV files here...
]

# DataFrame to combine all files
combined_df = pd.DataFrame()

# Process each parsed CSV file
for parsed_file in parsed_files:
    sample_id, group = add_sample_metadata(parsed_file)
    print(f"Processing {parsed_file} with Sample ID: {sample_id}, Group: {group}")

    # Load the parsed CSV file
    df = pd.read_csv(parsed_file)

    # Add Sample ID and Group columns
    df['Sample_ID'] = sample_id
    df['Group'] = group

    # Append to the combined DataFrame
    combined_df = pd.concat([combined_df, df], ignore_index=True)

# Save the final combined DataFrame to a new CSV file
combined_csv_path = "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/combined_final.csv"
combined_df.to_csv(combined_csv_path, index=False)
print(f"Final combined file saved as: {combined_csv_path}")
