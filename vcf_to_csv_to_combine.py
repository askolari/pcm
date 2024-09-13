import os
import pandas as pd
import cyvcf2
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

# Define directories
indels_dir = "/mnt/c/Users/agsko/dev/pcm/all_indels/"
snvs_dir = "/mnt/c/Users/agsko/dev/pcm/all_snvs/"

# Initialize an empty dataframe to hold all variants
all_variants = pd.DataFrame()

def is_vcf_empty(vcf_file):
    """
    Check if a VCF file contains any variant data by looking for lines after the header.
    """
    with open(vcf_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                return False  # Found data after headers, so it's not empty
    return True  # No data after headers, so it's empty

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
    if 'indel' in filename.lower():
        variant_type = 'indel'
    elif 'snv' in filename.lower():
        variant_type = 'snv'
    else:
        raise ValueError(f"Could not determine variant type (indel or snv) from filename: {filename}")
    
    return sample_id, sample_group, variant_type

def parse_vcf(vcf_file, sample_id, sample_group, variant_type):
    vcf = cyvcf2.VCF(vcf_file)
    variants = []
    
    for variant in vcf:
        info_dict = variant.INFO
        csq_field = info_dict.get('CSQ')
        if csq_field:
            csq_entries = csq_field.split(',')
            for csq in csq_entries:
                annotations = csq.split('|')
                filter_value = ','.join(variant.FILTER) if variant.FILTER is not None else "PASS"  # Handle non-iterable FILTER
                variant_data = {
                    'CHROM': variant.CHROM,
                    'POS': variant.POS,
                    'ID': variant.ID,
                    'REF': variant.REF,
                    'ALT': ','.join(variant.ALT),
                    'QUAL': variant.QUAL,
                    'FILTER': filter_value,
                    'INFO': variant.INFO,
                    'CSQ_Allele': annotations[0],
                    'CSQ_Consequence': annotations[1],
                    'CSQ_IMPACT': annotations[2],
                    'CSQ_SYMBOL': annotations[3],
                    'CSQ_Gene': annotations[4],
                    # Add more CSQ fields as needed
                    'Sample_ID': sample_id,
                    'Sample_Group': sample_group,
                    'Variant_Type': variant_type
                }
                variants.append(variant_data)
    
    return pd.DataFrame(variants)

# Iterate over files in the indels directory
for filename in os.listdir(indels_dir):
    vcf_file = os.path.join(indels_dir, filename)
    if os.path.isfile(vcf_file):  # Ensure it's a file, not a directory
        if is_vcf_empty(vcf_file):
            print(f"Warning: The VCF file {vcf_file} is empty. Skipping...")
            continue  # Skip this file since it contains no data
        try:
            sample_id, sample_group, variant_type = extract_sample_info_from_path(vcf_file)
        except ValueError as e:
            print(e)  # Log the error and skip this file
            continue
        variants_df = parse_vcf(vcf_file, sample_id, sample_group, variant_type)
        all_variants = pd.concat([all_variants, variants_df], ignore_index=True)

# Iterate over files in the snvs directory
for filename in os.listdir(snvs_dir):
    vcf_file = os.path.join(snvs_dir, filename)
    if os.path.isfile(vcf_file):  # Ensure it's a file, not a directory
        if is_vcf_empty(vcf_file):
            print(f"Warning: The VCF file {vcf_file} is empty. Skipping...")
            continue  # Skip this file since it contains no data
        try:
            sample_id, sample_group, variant_type = extract_sample_info_from_path(vcf_file)
        except ValueError as e:
            print(e)  # Log the error and skip this file
            continue
        variants_df = parse_vcf(vcf_file, sample_id, sample_group, variant_type)
        all_variants = pd.concat([all_variants, variants_df], ignore_index=True)

# Save the combined data to a CSV file with sample ID, group, and variant type information
all_variants.to_csv("combined_variants_with_sample_info.csv", index=False)
