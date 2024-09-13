#!/bin/bash

# Define the list of samples
samples=("Exome-seq_Control_1_S1" "Exome-seq_Control_2_S2" "Exome-seq_Sample_1_S3" "Exome-seq_Sample_2_S4" "Exome-seq_Sample_3_S5" "Exome-seq_Sample_4_S6" "Exome-seq_Sample_5_S7" "Exome-seq_Sample_6_S8" "Exome-seq_Sample_7_S9" "Exome-seq_Sample_8_S10" "Exome-seq_Sample_9_S11" "Exome-seq_Sample_10_S12" "Exome-seq_Sample_11_S13" "Exome-seq_Sample_12_S14" "Exome-seq_Sample_13_S15")

# Directory paths
vcf_dir="/mnt/c/Users/agsko/dev/pcm/all_indels"
output_dir="/mnt/c/Users/agsko/dev/pcm/all_indels"

# Loop over each sample
for sample in "${samples[@]}"; do
    # Generate Class I mutant peptides (shorter peptides, 17 flanking amino acids)
    pvacseq generate_protein_fasta "${vcf_dir}/${sample}_indels_filtered_lifted_annotated.vep.vcf" 17 "${output_dir}/${sample}_class1_mutant_peptides.txt" --sample-name DUMMY --mutant-only --pass-only
    
    # Generate Class II mutant peptides (longer peptides, 29 flanking amino acids)
    pvacseq generate_protein_fasta "${vcf_dir}/${sample}_indels_filtered_lifted_annotated.vep.vcf" 29 "${output_dir}/${sample}_class2_mutant_peptides.txt" --sample-name DUMMY --mutant-only --pass-only
done
