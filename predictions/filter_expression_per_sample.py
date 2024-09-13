
import pandas as pd

# Load the file
file_path = "/mnt/c/Users/agsko/dev/pcm/predictions/filtered_genes_with_counts_and_conditions.csv"
df = pd.read_csv(file_path)

# Define sample mappings
sample_mapping = {
    "X88": "S4",
    "X89": "S5",
    "X90": "S6",
    "X92": "S7",
    "X93": "S8",
    "X94": "S9",
    "X97": "S11",
    "X99": "S13",
    "X101": "S14"
}

# Loop through each sample ID and create a new file with the added RNA_Sample_ID column
for sample_id, rna_sample_id in sample_mapping.items():
    # Filter the dataframe for the current sample
    sample_df = df[df['Sample.ID'] == sample_id]
    
    # Add the RNA_Sample_ID column
    sample_df['Sample_ID'] = rna_sample_id
    
    # Save the filtered dataframe to a new file
    output_file = f"/mnt/c/Users/agsko/dev/pcm/predictions/{sample_id}_filtered_counts.csv"
    sample_df.to_csv(output_file, index=False)

    print(f"File saved for {sample_id} as {output_file}")

import pandas as pd
import os

# Define file paths for matched peptides and filtered counts
matched_class2_file = "/mnt/c/Users/agsko/dev/pcm/predictions/matched_class2_all_peptides.csv"
matched_class1_file = "/mnt/c/Users/agsko/dev/pcm/predictions/matched_class1_all_peptides.csv"
filtered_counts_files = [
    "/mnt/c/Users/agsko/dev/pcm/predictions/X93_filtered_counts.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X94_filtered_counts.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X97_filtered_counts.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X99_filtered_counts.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X101_filtered_counts.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X88_filtered_counts.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X89_filtered_counts.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X90_filtered_counts.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X92_filtered_counts.csv"
]

# Load the matched peptides files
matched_class2_df = pd.read_csv(matched_class2_file)
matched_class1_df = pd.read_csv(matched_class1_file)

# Loop through filtered counts files and check for matches
for file_path in filtered_counts_files:
    filtered_counts_df = pd.read_csv(file_path)
    
    # Perform matching based on 'Gene.name' and 'Sample_ID'
    matched_class2_filtered = pd.merge(matched_class2_df, filtered_counts_df, left_on=['Gene', 'Sample_ID'], right_on=['Gene.name', 'Sample_ID'], how='inner')
    matched_class1_filtered = pd.merge(matched_class1_df, filtered_counts_df, left_on=['Gene', 'Sample_ID'], right_on=['Gene.name', 'Sample_ID'], how='inner')
    
    # Save the matches to new files
    if not matched_class2_filtered.empty:
        output_file_class2 = f"/mnt/c/Users/agsko/dev/pcm/predictions/{os.path.basename(file_path).replace('_filtered_counts.csv', '')}_matched_class2.csv"
        matched_class2_filtered.to_csv(output_file_class2, index=False)
        print(f"Matched rows saved to {output_file_class2}")
    
    if not matched_class1_filtered.empty:
        output_file_class1 = f"/mnt/c/Users/agsko/dev/pcm/predictions/{os.path.basename(file_path).replace('_filtered_counts.csv', '')}_matched_class1.csv"
        matched_class1_filtered.to_csv(output_file_class1, index=False)
        print(f"Matched rows saved to {output_file_class1}")

