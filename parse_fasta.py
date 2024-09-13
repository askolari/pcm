import pandas as pd
import os

# List of file paths
file_paths = [
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_12_S14_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_1_S3_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_1_S3_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_1_S3_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_1_S3_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_2_S4_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_2_S4_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_2_S4_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_2_S4_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_3_S5_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_3_S5_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_3_S5_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_3_S5_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_4_S6_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_4_S6_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_4_S6_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_4_S6_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_5_S7_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_5_S7_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_5_S7_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_5_S7_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_6_S8_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_6_S8_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_6_S8_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_6_S8_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_7_S9_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_7_S9_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_7_S9_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_7_S9_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_8_S10_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_8_S10_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_8_S10_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_8_S10_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_9_S11_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_9_S11_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_9_S11_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_9_S11_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_10_S12_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_10_S12_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_10_S12_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_10_S12_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_11_S13_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_11_S13_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_11_S13_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_11_S13_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_13_S15_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_13_S15_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_13_S15_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_13_S15_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_12_S14_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_12_S14_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.fasta",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_12_S14_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.fasta"
]

# Create an empty list to store dataframes
dfs = []

# Loop through each file and process
for file_path in file_paths:
    # Load the file without headers
    df = pd.read_csv(file_path, header=None)
    
    # Check the number of columns in the DataFrame
    if df.shape[1] == 1:
        # If only one column, treat it as raw information and sequence in alternating rows
        df['Raw_Info'] = df[0][::2].values  # Take every other row for Raw_Info
        df['Sequence'] = df[0][1::2].values  # Take alternating rows for Sequence
        df = df.drop(columns=[0])  # Drop the original single column
    else:
        # If two columns already exist, just assign column names
        df.columns = ['Raw_Info', 'Sequence']
    
    # Split the 'Raw_Info' column into multiple columns
    df[['MT', 'Index', 'Gene', 'Transcript', 'Transcript_Version', 'Type', 'Position']] = df['Raw_Info'].str.split('.', expand=True)
    
    # Add a column for the sample ID extracted from the file name
    sample_id = os.path.basename(file_path).split('_')[3]
    df['Sample_ID'] = sample_id
    
    # Reorder columns for better readability
    df = df[['MT', 'Index', 'Gene', 'Transcript', 'Transcript_Version', 'Type', 'Position', 'Sequence', 'Sample_ID']]
    
    # Append to the list of dataframes
    dfs.append(df)

# Concatenate all the dataframes into one
combined_df = pd.concat(dfs, ignore_index=True)

# Save the combined DataFrame to a new CSV file
output_file_path = '/mnt/c/Users/agsko/dev/pcm/predictions/combined_final_peptides.csv'
combined_df.to_csv(output_file_path, index=False)

print(f"Combined dataset saved to {output_file_path}")