import pandas as pd
from collections import defaultdict

# Define file paths
file_paths = [
    "/mnt/c/Users/agsko/dev/pcm/predictions/X101_matched_class1.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X88_matched_class1.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X89_matched_class1.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X90_matched_class1.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X92_matched_class1.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X93_matched_class1.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X94_matched_class1.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X97_matched_class1.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X99_matched_class1.csv"
]

# Dictionary to hold peptides and the Sample_IDs they appear in across different samples
peptide_samples = defaultdict(set)  # Store as set to avoid duplicate Sample_IDs

# List to hold the merged data for output
merged_data = []

# Read each file, extract the 'peptide' column, and keep track of where each peptide is found with different Sample_IDs
for file_path in file_paths:
    sample_name = file_path.split('/')[-1].replace('_matched_class1.csv', '')  # Get the sample name from the file path
    try:
        df = pd.read_csv(file_path)
        for _, row in df.iterrows():
            peptide = row['peptide']
            sample_id = row['Sample_ID']
            peptide_samples[peptide].add(sample_id)  # Add unique Sample_ID for each peptide
            merged_data.append({**row.to_dict(), 'Sample_Found': sample_name})  # Add the row data with the sample name
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        continue

# Filter peptides that appear in at least two samples with different Sample_IDs
common_peptides = {peptide: sample_ids for peptide, sample_ids in peptide_samples.items() if len(sample_ids) >= 2}

# Create a dataframe for the common peptides
common_peptides_data = [row for row in merged_data if row['peptide'] in common_peptides]

# Convert to a DataFrame and add the 'Samples_Found_In' column
output_df = pd.DataFrame(common_peptides_data)
output_df['Samples_Found_In'] = output_df['peptide'].apply(lambda p: ', '.join(common_peptides[p]))

# Save to a new file
output_file = "/mnt/c/Users/agsko/dev/pcm/predictions/common_peptides_in_different_samples.csv"
output_df.to_csv(output_file, index=False)

print(f"Common peptides saved to {output_file}")
