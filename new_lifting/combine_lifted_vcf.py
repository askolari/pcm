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

# Directory containing the lifted VCF files
vcf_dir = "/mnt/c/Users/agsko/dev/pcm/new_lifting"

# Initialize an empty dataframe to hold all variants
all_variants = []

def extract_sample_info_from_path(vcf_file):
    """
    Extract the sample ID, group, and variant type from the VCF file path.
    Assumes the sample ID is in the format '_S1_' within the filename.
    """
    filename = os.path.basename(vcf_file)
    
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

def parse_vcf(vcf_file, sample_id, sample_group, variant_type):
    """
    Parse the VCF file to extract the relevant columns and add sample metadata.
    """
    # Initialize lists to hold the parsed data
    records = []
    
    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = fields[1]
            id_ = fields[2]
            ref = fields[3]
            alt = fields[4]
            qual = fields[5]
            filter_ = fields[6]
            info = fields[7]
            csq_fields = info.split(';')[0].split('|')  # Assuming CSQ is the first field in INFO
            if len(csq_fields) > 4:  # Ensure there's enough CSQ data
                csq_allele = csq_fields[0]
                csq_consequence = csq_fields[1]
                csq_impact = csq_fields[2]
                csq_symbol = csq_fields[3]
                csq_gene = csq_fields[4]
            else:
                csq_allele = csq_consequence = csq_impact = csq_symbol = csq_gene = None
            records.append([
                chrom, pos, id_, ref, alt, qual, filter_, info,
                csq_allele, csq_consequence, csq_impact, csq_symbol, csq_gene,
                sample_id, sample_group, variant_type
            ])
    
    # Create a DataFrame for the current VCF file
    columns = [
        'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
        'CSQ_Allele', 'CSQ_Consequence', 'CSQ_IMPACT', 'CSQ_SYMBOL', 'CSQ_Gene',
        'Sample_ID', 'Sample_Group', 'Variant_Type'
    ]
    df = pd.DataFrame(records, columns=columns)
    return df

# Iterate over files in the VCF directory
for filename in os.listdir(vcf_dir):
    vcf_file = os.path.join(vcf_dir, filename)
    if os.path.isfile(vcf_file) and filename.endswith('.vcf'):  # Ensure it's a file and a VCF file
        try:
            sample_id, sample_group, variant_type = extract_sample_info_from_path(vcf_file)
        except ValueError as e:
            print(e)  # Log the error and skip this file
            continue
        variants_df = parse_vcf(vcf_file, sample_id, sample_group, variant_type)
        all_variants.append(variants_df)

# Combine all DataFrames into one
combined_variants = pd.concat(all_variants, ignore_index=True)

# Save the combined data to a CSV file with the specified columns
output_file = os.path.join(vcf_dir, "combined_lifted_vcf_variants.csv")
combined_variants.to_csv(output_file, index=False)

print(f"Combined VCF dataset saved to {output_file}")
