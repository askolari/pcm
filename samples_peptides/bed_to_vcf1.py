import pandas as pd

# Load the merged BED file
merged_bed_path = "/mnt/c/Users/agsko/dev/pcm/Exome-seq_Control_2_S2_snvs_filtered_lifted_with_info.bed"
output_vcf_path = "/mnt/c/Users/agsko/dev/pcm/Exome-seq_COntrol_2_S2_snvs_filtered_lifted.vcf"

merged_bed = pd.read_csv(merged_bed_path, sep='\t', header=None)
merged_bed.columns = ['CHROM', 'POS', 'END', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR']

# Create a new DataFrame for the VCF
vcf_df = merged_bed[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR']]

# Load the original VCF header
original_vcf_path = "/mnt/c/Users/agsko/dev/pcm/Exome-seq_Control_2_S2_snvs_filtered.vcf"
with open(original_vcf_path, 'r') as file:
    vcf_header = []
    for line in file:
        if line.startswith("#"):
            vcf_header.append(line.strip())
        else:
            break

# Save the VCF file
with open(output_vcf_path, 'w') as file:
    # Write the VCF header
    for line in vcf_header:
        file.write(line + '\n')
    # Write the VCF data
    vcf_df.to_csv(file, sep='\t', index=False, header=False)

print("VCF file has been created.")
