import pandas as pd

# Load the dataset with TPM values
final_variants_with_tpm_df = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/merged_files/final_variants_with_TPM.csv')

# Filter the rows where 'TPM.values' is not empty (NaN)
filtered_df = final_variants_with_tpm_df.dropna(subset=['TPM.values'])

# Save the filtered DataFrame to a new CSV file
filtered_df.to_csv('/mnt/c/Users/agsko/dev/pcm/merged_files/final_variants_with_TPM_filtered.csv', index=False)

print("Filtered dataset with non-empty TPM values has been saved.")

import pandas as pd

# Load the dataset with TPM values
final_variants_with_tpm_df = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/merged_files/final_variants_with_TPM.csv')

# Filter the rows where 'TPM.values' is not empty (NaN) and not equal to 0
filtered_df = final_variants_with_tpm_df.dropna(subset=['TPM.values'])
filtered_df = filtered_df[filtered_df['TPM.values'] != 0]

# Save the filtered DataFrame to a new CSV file
filtered_df.to_csv('/mnt/c/Users/agsko/dev/pcm/merged_files/final_variants_with_TPM_filtered_nonzero.csv', index=False)

print("Filtered dataset with non-zero TPM values has been saved.")

import pandas as pd

# Load the filtered dataset with non-zero TPM values
filtered_df = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/merged_files/final_variants_with_TPM_filtered_nonzero.csv')

# Get the unique Sample_IDs
sample_ids = filtered_df['Sample_ID'].unique()

# Create a separate dataset for each Sample_ID
for sample_id in sample_ids:
    sample_df = filtered_df[filtered_df['Sample_ID'] == sample_id]
    
    # Define the output file path
    output_file = f'/mnt/c/Users/agsko/dev/pcm/merged_files/final_variants_with_TPM_{sample_id}.csv'
    
    # Save the dataset for this Sample_ID
    sample_df.to_csv(output_file, index=False)
    
    print(f"Dataset for Sample_ID {sample_id} saved to {output_file}")

