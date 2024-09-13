import pandas as pd

# Paths to the BED files
original_bed_path = "/mnt/c/Users/agsko/dev/pcm/Exome-seq_Control_2_S2_snvs_filtered.bed"
lifted_simplified_bed_path = "/mnt/c/Users/agsko/dev/pcm/Exome-seq_Control_2_S2_snvs_filtered_lifted_simplified.bed"
output_merged_bed_path = "/mnt/c/Users/agsko/dev/pcm/Exome-seq_Control_2_S2_snvs_filtered_lifted_with_info.bed"

# Load the original BED file
original_bed = pd.read_csv(original_bed_path, sep='\t', header=None)
original_bed.columns = ['CHROM', 'START', 'END', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR']

# Load the lifted simplified BED file
lifted_bed = pd.read_csv(lifted_simplified_bed_path, sep='\t', header=None)
lifted_bed.columns = ['CHROM', 'START', 'END', 'ID']

# Adjust positions: BED files are 0-based while VCF files are 1-based, so adjust the positions accordingly
lifted_bed['START'] = lifted_bed['START'] + 1

# Ensure both dataframes have the same length
if len(original_bed) != len(lifted_bed):
    raise ValueError("The original and lifted BED files do not have the same number of entries.")

# Merge the coordinates from the lifted BED into the original BED
merged_bed = original_bed.copy()
merged_bed['CHROM'] = lifted_bed['CHROM']
merged_bed['START'] = lifted_bed['START']
merged_bed['END'] = lifted_bed['END']

# Save the merged BED file with all the variant information
merged_bed.to_csv(output_merged_bed_path, sep='\t', index=False, header=False)

print("Merged BED file has been created.")
