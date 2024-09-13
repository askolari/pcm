import pandas as pd

# Load the combined dataset
input_file = "/mnt/c/Users/agsko/dev/pcm/combined_all_samples_with_groups.csv"
output_file = "/mnt/c/Users/agsko/dev/pcm/all_tsv/final_variants_mm10.csv"

# Load the combined TSV data
df = pd.read_csv(input_file)

# Function to calculate ALT_Reads from the TUMOR.AU column
def calculate_alt_reads(tumor_au):
    if pd.isna(tumor_au):
        return 0
    try:
        return sum(map(int, tumor_au.split(',')))
    except Exception as e:
        print(f"Error processing TUMOR.AU value: {tumor_au}. Error: {e}")
        return 0

# Apply the function to calculate ALT_Reads
df['ALT_Reads'] = df['TUMOR.AU'].apply(calculate_alt_reads)

# Calculate VAF (Variant Allele Frequency)
df['VAF'] = df.apply(lambda row: row['ALT_Reads'] / row['TUMOR.DP'] if row['TUMOR.DP'] > 0 else 0, axis=1)

# Save the updated DataFrame with the VAF column
df.to_csv(output_file, index=False)

print(f"Combined dataset with VAF saved to {output_file}")
