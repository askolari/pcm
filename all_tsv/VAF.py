import os
import pandas as pd
import re

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

# Define the directory containing the TSV files
tsv_dir = "/mnt/c/Users/agsko/dev/pcm/all_tsv"


# Initialize an empty dataframe to hold all variants
all_variants = []

def is_tsv_empty(tsv_file):
    """
    Check if a TSV file contains any variant data by checking if it's empty.
    """
    df = pd.read_csv(tsv_file, sep='\t')
    return df.empty

def extract_sample_info_from_path(tsv_file):
    """
    Extract the sample ID, group, and variant type from the TSV file path.
    Assumes the sample ID is in the format '_S1_' within the filename.
    """
    filename = os.path.basename(tsv_file)
    
    # Extract sample ID using regex
    sample_id_match = re.search(r'_S\d{1,2}_', filename)  # Match patterns like '_S1_', '_S12_', etc.
    if sample_id_match:
        sample_id = sample_id_match.group(0).strip('_')  # Extract and clean the sample ID
    else:
        raise ValueError(f"Could not extract sample ID from filename: {filename}")
    
    # Get sample group from metadata
    sample_group = sample_metadata.get(sample_id)
    if not sample_group:
        raise ValueError(f"No group found for sample ID: {sample_id}")
    
    # Determine variant type (indel or snv)
    if 'indels' in filename.lower():
        variant_type = 'indel'
    elif 'snvs' in filename.lower():
        variant_type = 'snv'
    else:
        raise ValueError(f"Could not determine variant type (indel or snv) from filename: {filename}")
    
    return sample_id, sample_group, variant_type

def parse_tsv(tsv_file, sample_id, sample_group, variant_type):
    # Load TSV data
    df = pd.read_csv(tsv_file, sep='\t')

    # Add Sample_ID, Sample_Group, and Variant_Type to the dataframe
    df['Sample_ID'] = sample_id
    df['Sample_Group'] = sample_group
    df['Variant_Type'] = variant_type

    # Calculate VAF (Variant Allele Frequency) if columns for TUMOR.AU, TUMOR.DP, etc. exist
    if 'TUMOR.AU' in df.columns and 'TUMOR.DP' in df.columns:
        df['ALT_Reads'] = df['TUMOR.AU'].apply(lambda x: sum(map(int, x.split(','))))
        df['VAF'] = df['ALT_Reads'] / df['TUMOR.DP']
    
    return df

# Iterate over files in the TSV directory
for filename in os.listdir(tsv_dir):
    tsv_file = os.path.join(tsv_dir, filename)
    if os.path.isfile(tsv_file):  # Ensure it's a file, not a directory
        if is_tsv_empty(tsv_file):
            print(f"Warning: The TSV file {tsv_file} is empty. Skipping...")
            continue  # Skip this file since it contains no data
        try:
            sample_id, sample_group, variant_type = extract_sample_info_from_path(tsv_file)
        except ValueError as e:
            print(e)  # Log the error and skip this file
            continue
        variants_df = parse_tsv(tsv_file, sample_id, sample_group, variant_type)
        all_variants.append(variants_df)

# Combine all DataFrames into one
combined_variants = pd.concat(all_variants, ignore_index=True)

# Save the combined data to a CSV file with sample ID, group, and variant type information
output_file = os.path.join(tsv_dir, "combined_tsv_variants.csv")
combined_variants.to_csv(output_file, index=False)

print(f"Combined dataset with VAF saved to {output_file}")