import pandas as pd
import os

# Define directories
bed_dir = "/mnt/c/Users/agsko/dev/pcm/all_tsv/"
tsv_dir = "/mnt/c/Users/agsko/dev/pcm/all_tsv/"
output_dir = "/mnt/c/Users/agsko/dev/pcm/all_tsv/"

# List of file prefixes (these should match the prefixes of your files)
file_prefixes = [
    "Exome-seq_Sample_10_S12_indels_filtered",
    "Exome-seq_Sample_10_S12_snvs_filtered",
    "Exome-seq_Sample_11_S13_indels_filtered",
    "Exome-seq_Sample_11_S13_snvs_filtered",
    "Exome-seq_Sample_12_S14_indels_filtered",
    "Exome-seq_Sample_12_S14_snvs_filtered",
    "Exome-seq_Sample_13_S15_indels_filtered",
    "Exome-seq_Sample_13_S15_snvs_filtered",
    "Exome-seq_Control_1_S1_indels_filtered",
    "Exome-seq_Control_1_S1_snvs_filtered",
    "Exome-seq_Control_2_S2_indels_filtered",
    "Exome-seq_Control_2_S2_snvs_filtered",
    "Exome-seq_Sample_1_S3_indels_filtered",
    "Exome-seq_Sample_1_S3_snvs_filtered",
    "Exome-seq_Sample_2_S4_indels_filtered",
    "Exome-seq_Sample_2_S4_snvs_filtered",
    "Exome-seq_Sample_3_S5_indels_filtered",
    "Exome-seq_Sample_3_S5_snvs_filtered",
    "Exome-seq_Sample_4_S6_indels_filtered",
    "Exome-seq_Sample_4_S6_snvs_filtered",
    "Exome-seq_Sample_5_S7_indels_filtered",
    "Exome-seq_Sample_5_S7_snvs_filtered",
    "Exome-seq_Sample_6_S8_indels_filtered",
    "Exome-seq_Sample_6_S8_snvs_filtered",
    "Exome-seq_Sample_7_S9_indels_filtered",
    "Exome-seq_Sample_7_S9_snvs_filtered",
    "Exome-seq_Sample_8_S10_indels_filtered",
    "Exome-seq_Sample_8_S10_snvs_filtered",
    "Exome-seq_Sample_9_S11_indels_filtered",
    "Exome-seq_Sample_9_S11_snvs_filtered"
]

for file_prefix in file_prefixes:
    try:
        # Construct full paths based on the correct filenames
        lifted_simplified_bed_path = os.path.join(bed_dir, f"{file_prefix}_simplified_tsv_lifted.bed")
        original_tsv_path = os.path.join(tsv_dir, f"{file_prefix}.tsv")
        output_tsv_path = os.path.join(output_dir, f"{file_prefix}_lifted.tsv")

        # Load the simplified lifted BED file
        lifted_bed = pd.read_csv(lifted_simplified_bed_path, sep='\t', header=None)
        lifted_bed.columns = ['CHROM', 'POS', 'END', 'ID']

        # Adjust positions: BED files are 0-based while VCF/TSV files are 1-based, so adjust the positions accordingly
        lifted_bed['POS'] = lifted_bed['POS'] + 1

        # Load the original TSV data
        original_tsv = pd.read_csv(original_tsv_path, sep='\t')
        
        # Ensure both dataframes have the same length
        if len(original_tsv) != len(lifted_bed):
            raise ValueError(f"The original TSV and lifted BED files do not have the same number of entries for {file_prefix}.")

        # Merge the coordinates from the lifted BED into the original TSV
        merged_tsv = original_tsv.copy()
        merged_tsv['CHROM'] = lifted_bed['CHROM']
        merged_tsv['POS'] = lifted_bed['POS']

        # Save the TSV file with the updated coordinates
        merged_tsv.to_csv(output_tsv_path, sep='\t', index=False)

        print(f"Lifted TSV file created: {output_tsv_path}")
    except Exception as e:
        print(f"An error occurred while processing {file_prefix}: {e}")
