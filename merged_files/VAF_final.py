import pandas as pd

# Load the CSV file into a pandas DataFrame
df = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/merged_files/combined_all_samples_with_groups.csv')

# Convert relevant columns to numeric, forcing errors to NaN
numeric_columns = ['TUMOR.DP', 'TUMOR.AU', 'TUMOR.CU', 'TUMOR.GU', 'TUMOR.TU']
df[numeric_columns] = df[numeric_columns].apply(pd.to_numeric, errors='coerce')

def calculate_vaf(row):
    # Identify the correct column for the alternate allele
    if row['ALT_x'] == 'A':
        alt_reads = row['TUMOR.AU']
    elif row['ALT_x'] == 'C':
        alt_reads = row['TUMOR.CU']
    elif row['ALT_x'] == 'G':
        alt_reads = row['TUMOR.GU']
    elif row['ALT_x'] == 'T':
        alt_reads = row['TUMOR.TU']
    else:
        alt_reads = None  # Handle unexpected cases

    # Print the current row's ALT_x value and the corresponding alt_reads value
    print(f"ALT_x: {row['ALT_x']}, alt_reads: {alt_reads}")

    # Check if alt_reads is NaN or if TUMOR.DP is 0 or NaN
    if pd.isna(alt_reads) or row['TUMOR.DP'] == 0 or pd.isna(row['TUMOR.DP']):
        return None
    
    # Calculate VAF
    return alt_reads / row['TUMOR.DP']

# Apply the function to calculate VAF for each row
df['VAF'] = df.apply(calculate_vaf, axis=1)

# Inspect the first few rows to ensure VAF is calculated correctly
print(df[['CHROM_x', 'POS', 'REF_x', 'ALT_x', 'TUMOR.DP', 'VAF']].head())

# Save the updated DataFrame with the VAF column to a new CSV file
output_file_path = '/mnt/c/Users/agsko/dev/pcm/merged_files/combined_all_samples_with_groups_with_VAF.csv'
df.to_csv(output_file_path, index=False)

print(f"File saved successfully to {output_file_path}")
