import pandas as pd

# Define file paths
combined_file_path = '/mnt/c/Users/agsko/dev/pcm/combined_all_samples.csv'
output_file_path = '/mnt/c/Users/agsko/dev/pcm/combined_all_samples_with_groups.csv'

# Load the combined dataset
combined_df = pd.read_csv(combined_file_path)

# Convert Sample_ID to string and prepend 'S'
combined_df['Sample_ID'] = 'S' + combined_df['Sample_ID'].astype(str)

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

# Assign Sample_Group based on Sample_ID
combined_df['Sample_Group'] = combined_df['Sample_ID'].map(sample_metadata)

# Verify if any Sample_IDs didn't match
unmatched_samples = combined_df[combined_df['Sample_Group'].isnull()]['Sample_ID'].unique()
if len(unmatched_samples) > 0:
    print(f"Unmatched Sample_IDs: {unmatched_samples}")
else:
    print("All Sample_IDs successfully mapped to Sample_Groups.")

# Save the updated DataFrame
combined_df.to_csv(output_file_path, index=False)

print(f"Updated combined file with Sample_Group saved as: {output_file_path}")

