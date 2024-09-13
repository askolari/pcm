import os
import pandas as pd

# Function to parse the VCF file
def parse_vcf(vcf_file):
    data = []
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
            data.append([chrom, pos, id_, ref, alt, qual, filter_, info])
    
    columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    df = pd.DataFrame(data, columns=columns)
    return df

# Function to add sample information
def add_sample_info(df, sample_id, sample_group):
    df['Sample_ID'] = sample_id
    df['Sample_Group'] = sample_group
    return df

# Function to extract sample ID from the file path
def extract_sample_id(file_path):
    filename = os.path.basename(file_path)
    return filename.split('_')[3]  # Extracts the Sample ID (e.g., S1, S2, etc.)

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

# List of VCF file paths
vcf_files = [
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_12_S14_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_12_S14_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_13_S15_indel_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_13_S15_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Control_1_S1_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Control_2_S2_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_1_S3_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_1_S3_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_2_S4_indel_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_2_S4_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_3_S5_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_3_S5_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_4_S6_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_4_S6_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_5_S7_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_5_S7_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_6_S8_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_6_S8_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_7_S9_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_7_S9_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_8_S10_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_8_S10_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_9_S11_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_9_S11_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_10_S12_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_10_S12_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_11_S13_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_11_S13_snvs_filtered_lifted_annotated.vep.vcf",
]



combined_df = pd.DataFrame()

for vcf_file in vcf_files:
    # Parse the VCF file
    df = parse_vcf(vcf_file)
    
    # Extract sample ID from the filename
    sample_id = extract_sample_id(vcf_file)
    
    # Get the sample group from the metadata
    sample_group = sample_metadata.get(sample_id, 'Unknown')
    
    # Add sample information to the DataFrame
    df = add_sample_info(df, sample_id, sample_group)
    
    # Append the DataFrame to the combined DataFrame
    combined_df = pd.concat([combined_df, df], ignore_index=True)

# Save the combined DataFrame to a CSV file
output_file = '/mnt/c/Users/agsko/dev/pcm/combined_alterations_with_metadata.csv'
combined_df.to_csv(output_file, index=False)

print(f"Combined CSV saved as: {output_file}")
