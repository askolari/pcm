import pandas as pd

# Load the original BED file
original_bed = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/Exome-seq_Control_2_S2_snvs_filtered.bed', sep='\t', header=None)

# Create a simplified BED file with only the necessary columns for liftover
simplified_bed = original_bed.iloc[:, [0, 1, 2, 3]]
simplified_bed.columns = ['CHROM', 'START', 'END', 'ID']

# Save the simplified BED file
simplified_bed.to_csv('/mnt/c/Users/agsko/dev/pcm/Exome-seq_Control_2_S2_snvs_filtered_simplified.bed', sep='\t', index=False, header=False)
