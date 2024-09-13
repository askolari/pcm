import pandas as pd
import os

# Define directories
bed_dir = "/mnt/c/Users/agsko/dev/pcm/new_lifting/"
vcf_dir = "/mnt/c/Users/agsko/dev/pcm/new_lifting/"
output_dir = "/mnt/c/Users/agsko/dev/pcm/new_lifting/"

# List of file prefixes (these should match the prefixes of your files)
file_prefixes = [
    "Exome-seq_Sample_13_S15_snvs_filtered",
    "Exome-seq_Sample_13_S15_indels_filtered",
    "Exome-seq_Sample_12_S14_snvs_filtered",
    "Exome-seq_Sample_12_S14_indels_filtered",
    "Exome-seq_Sample_11_S13_snvs_filtered",
    "Exome-seq_Sample_11_S13_indels_filtered",
    "Exome-seq_Sample_10_S12_snvs_filtered",
    "Exome-seq_Sample_10_S12_indels_filtered",
    "Exome-seq_Sample_9_S11_snvs_filtered",
    "Exome-seq_Sample_9_S11_indels_filtered",
    "Exome-seq_Sample_8_S10_snvs_filtered",
    "Exome-seq_Sample_8_S10_indels_filtered",
    "Exome-seq_Sample_7_S9_snvs_filtered",
    "Exome-seq_Sample_7_S9_indels_filtered",
    "Exome-seq_Sample_6_S8_snvs_filtered",
    "Exome-seq_Sample_6_S8_indels_filtered",
    "Exome-seq_Sample_5_S7_snvs_filtered",
    "Exome-seq_Sample_5_S7_indels_filtered",
    "Exome-seq_Sample_4_S6_snvs_filtered",
    "Exome-seq_Sample_4_S6_indels_filtered",
    "Exome-seq_Sample_3_S5_snvs_filtered",
    "Exome-seq_Sample_3_S5_indels_filtered",
    "Exome-seq_Sample_2_S4_snvs_filtered",
    "Exome-seq_Sample_2_S4_indels_filtered",
    "Exome-seq_Sample_1_S3_snvs_filtered",
    "Exome-seq_Sample_1_S3_indels_filtered",
    "Exome-seq_Control_2_S2_snvs_filtered",
    "Exome-seq_Control_2_S2_indels_filtered",
    "Exome-seq_Control_1_S1_snvs_filtered",
    "Exome-seq_Control_1_S1_indels_filtered"
]

for file_prefix in file_prefixes:
    try:
        # Construct full paths based on the correct filenames
        lifted_simplified_bed_path = os.path.join(bed_dir, f"{file_prefix}_simplified_vcf_lifted.bed")
        original_vcf_path = os.path.join(vcf_dir, f"{file_prefix}.vcf")
        output_vcf_path = os.path.join(output_dir, f"{file_prefix}_lifted.vcf")

        # Load the simplified lifted BED file
        lifted_bed = pd.read_csv(lifted_simplified_bed_path, sep='\t', header=None)
        lifted_bed.columns = ['CHROM', 'POS', 'END', 'ID']

        # Adjust positions: BED files are 0-based while VCF files are 1-based, so adjust the positions accordingly
        lifted_bed['POS'] = lifted_bed['POS'] + 1

        # Load the original VCF data (excluding the header)
        original_vcf = pd.read_csv(original_vcf_path, sep='\t', comment='#', header=None)
        original_vcf.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR']

        # Ensure both dataframes have the same length
        if len(original_vcf) != len(lifted_bed):
            raise ValueError(f"The original VCF and lifted BED files do not have the same number of entries for {file_prefix}.")

        # Merge the coordinates from the lifted BED into the original VCF
        merged_vcf = original_vcf.copy()
        merged_vcf['CHROM'] = lifted_bed['CHROM']
        merged_vcf['POS'] = lifted_bed['POS']

        # Load the original VCF header
        with open(original_vcf_path, 'r') as file:
            vcf_header = []
            for line in file:
                if line.startswith("#"):
                    vcf_header.append(line.strip())
                else:
                    break

        # Save the VCF file with the updated coordinates
        with open(output_vcf_path, 'w') as file:
            # Write the VCF header
            for line in vcf_header:
                file.write(line + '\n')
            # Write the updated VCF data
            merged_vcf.to_csv(file, sep='\t', index=False, header=False)

        print(f"Lifted VCF file created: {output_vcf_path}")
    except Exception as e:
        print(f"An error occurred while processing {file_prefix}: {e}")



