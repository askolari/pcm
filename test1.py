import pandas as pd
import os

# Directory to save the parsed output files
output_dir = "/mnt/c/Users/agsko/dev/pcm/predictions"

# Function to parse the 'Header' column and split it by periods
def parse_header_column(file_path, output_dir):
    # Load the file
    df = pd.read_csv(file_path, delimiter="\t", header=None, names=["Header", "Sequence"])
    
    # Split the 'Header' column by '.' and expand into new columns
    header_split = df['Header'].str.split('.', expand=True)

    # Dynamically generate column names for the split parts
    header_split.columns = [f'Column{i+1}' for i in range(header_split.shape[1])]

    # Concatenate the original dataframe with the new columns
    df_parsed = pd.concat([df, header_split], axis=1)

    # Define the output file name and save the parsed DataFrame
    output_file_name = os.path.basename(file_path).replace(".fasta", "_parsed.csv")
    output_path = os.path.join(output_dir, output_file_name)
    df_parsed.to_csv(output_path, index=False)
    print(f"Parsed and saved: {output_path}")

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



# Loop through each file and apply the function
for file_path in file_paths:
    parse_header_column(file_path, output_dir)
