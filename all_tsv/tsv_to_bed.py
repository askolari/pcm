import os
import pandas as pd

def tsv_to_bed_simplified(tsv_file, bed_file):
    df = pd.read_csv(tsv_file, sep='\t')

    # Ensure the necessary columns exist in the dataframe
    if 'CHROM' not in df.columns or 'POS' not in df.columns:
        raise ValueError("TSV file is missing one or more required columns: 'CHROM', 'POS'")

    # Create a simplified BED file with only the necessary columns
    simplified_bed = df[['CHROM', 'POS']].copy()
    simplified_bed['START'] = simplified_bed['POS'] - 1
    simplified_bed['END'] = simplified_bed['POS']
    simplified_bed['ID'] = df['REF']  # Assuming 'REF' to be the identifier

    # Rearrange columns to match the BED format
    simplified_bed = simplified_bed[['CHROM', 'START', 'END', 'ID']]

    # Save to a BED file
    simplified_bed.to_csv(bed_file, sep='\t', index=False, header=False)

# Define directory and file paths
tsv_dir = '/mnt/c/Users/agsko/dev/pcm/all_tsv/'
tsv_files = [
    "Exome-seq_Sample_1_S3_indels_filtered.tsv",
    "Exome-seq_Sample_12_S14_snvs_filtered.tsv",
    "Exome-seq_Sample_13_S15_indels_filtered.tsv",
    "Exome-seq_Sample_13_S15_snvs_filtered.tsv",
    "Exome-seq_Sample_10_S12_snvs_filtered.tsv",
    "Exome-seq_Sample_11_S13_indels_filtered.tsv",
    "Exome-seq_Sample_11_S13_snvs_filtered.tsv",
    "Exome-seq_Sample_12_S14_indels_filtered.tsv",
    "Exome-seq_Sample_8_S10_indels_filtered.tsv",
    "Exome-seq_Sample_8_S10_snvs_filtered.tsv",
    "Exome-seq_Sample_9_S11_indels_filtered.tsv",
    "Exome-seq_Sample_9_S11_snvs_filtered.tsv",
    "Exome-seq_Sample_10_S12_indels_filtered.tsv",
    "Exome-seq_Sample_6_S8_indels_filtered.tsv",
    "Exome-seq_Sample_6_S8_snvs_filtered.tsv",
    "Exome-seq_Sample_7_S9_indels_filtered.tsv",
    "Exome-seq_Sample_7_S9_snvs_filtered.tsv",
    "Exome-seq_Sample_3_S5_snvs_filtered.tsv",
    "Exome-seq_Sample_4_S6_indels_filtered.tsv",
    "Exome-seq_Sample_4_S6_snvs_filtered.tsv",
    "Exome-seq_Sample_5_S7_indels_filtered.tsv",
    "Exome-seq_Sample_5_S7_snvs_filtered.tsv",
    "Exome-seq_Sample_1_S3_snvs_filtered.tsv",
    "Exome-seq_Sample_2_S4_indels_filtered.tsv",
    "Exome-seq_Sample_2_S4_snvs_filtered.tsv",
    "Exome-seq_Sample_3_S5_indels_filtered.tsv",
    "Exome-seq_Control_1_S1_indels_filtered.tsv",
    "Exome-seq_Control_1_S1_snvs_filtered.tsv",
    "Exome-seq_Control_2_S2_indels_filtered.tsv",
    "Exome-seq_Control_2_S2_snvs_filtered.tsv"
]

# Process each TSV file
for tsv_file in tsv_files:
    tsv_path = os.path.join(tsv_dir, tsv_file)
    bed_file = tsv_path.replace('.tsv', '_simplified.bed')
    tsv_to_bed_simplified(tsv_path, bed_file)
    print(f"Converted {tsv_file} to BED format.")
