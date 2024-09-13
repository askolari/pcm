import pandas as pd
import os

# List of sample IDs
sample_ids = ['S3','S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15']

# Base file path for the predictions
base_file_path = '/mnt/c/Users/agsko/dev/pcm/predictions/'

# Loop through each sample ID and process the corresponding file
for sample_id in sample_ids:
    # Construct the file path
    input_file = os.path.join(base_file_path, f'{sample_id}_class2_predictions.csv')
    
    # Load the peptide prediction dataset
    predictions_df = pd.read_csv(input_file)
    
    # Primary filter: IC50 < 500 nM
    high_affinity_peptides = predictions_df[predictions_df['ic50'] < 500]
    
    # Secondary filter: Percentile rank < 1%
    high_affinity_peptides = high_affinity_peptides[high_affinity_peptides['rank'] < 1]
    
    # Sort by IC50 value and Percentile rank (for better prioritization)
    prioritized_peptides = high_affinity_peptides.sort_values(by=['ic50', 'rank'])
    
    # Save the prioritized peptides to a new CSV file
    output_file = os.path.join(base_file_path, f'prioritized_peptides_{sample_id}.csv')
    prioritized_peptides.to_csv(output_file, index=False)
    
    print(f"Prioritized peptides saved to {output_file}")

import pandas as pd
import os

# List of sample IDs
sample_ids = ['S3','S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15']

# Base file path for the predictions
base_file_path = '/mnt/c/Users/agsko/dev/pcm/predictions/'

# Create an empty list to hold the dataframes for each sample
high_affinity_dfs = []

# Loop through each sample ID and process the corresponding file
for sample_id in sample_ids:
    # Construct the file path
    input_file = os.path.join(base_file_path, f'{sample_id}_class2_predictions.csv')
    
    # Load the peptide prediction dataset
    predictions_df = pd.read_csv(input_file)
    
    # Filter: IC50 < 50 nM
    high_affinity_peptides = predictions_df[predictions_df['ic50'] < 50]

    # Secondary filter: Percentile rank < 1%
    high_affinity_peptides = high_affinity_peptides[high_affinity_peptides['rank'] < 1]
    
    # Add a column to track the sample ID
    high_affinity_peptides['Sample_ID'] = sample_id
    
    # Save the high-affinity peptides for this sample to a new CSV file
    output_file = os.path.join(base_file_path, f'high_affinity_peptides_{sample_id}.csv')
    high_affinity_peptides.to_csv(output_file, index=False)
    
    # Append the dataframe to the list
    high_affinity_dfs.append(high_affinity_peptides)
    
    print(f"High-affinity peptides for {sample_id} saved to {output_file}")

# Combine all high-affinity peptides into a single dataframe
all_high_affinity_peptides = pd.concat(high_affinity_dfs, ignore_index=True)

# Save the combined dataframe to a new CSV file
combined_output_file = os.path.join(base_file_path, 'high_peptides_class2_combined_all_samples.csv')
all_high_affinity_peptides.to_csv(combined_output_file, index=False)

print(f"Combined high-affinity peptides saved to {combined_output_file}")
