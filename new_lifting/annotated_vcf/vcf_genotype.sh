#!/bin/bash

# Define the base directory where your VCF files are located
BASE_DIR="/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf"

# Define an array of the VCF filenames you want to process
vcf_files=(
    "Exome-seq_Sample_1_S3_snvs_filtered_lifted_annotated_corrected.vep.vcf"
    "Exome-seq_Sample_1_S3_indels_filtered_lifted_annotated_corrected.vep.vcf"
    "Exome-seq_Sample_2_S4_snvs_filtered_lifted_annotated_corrected.vep.vcf"
    "Exome-seq_Sample_2_S4_indels_filtered_lifted_annotated_corrected.vep.vcf"
    "Exome-seq_Sample_3_S5_snvs_filtered_lifted_annotated_corrected.vep.vcf"
    "Exome-seq_Sample_3_S5_indels_filtered_lifted_annotated_corrected.vep.vcf"
    "Exome-seq_Sample_4_S6_snvs_filtered_lifted_annotated_corrected.vep.vcf"
    "Exome-seq_Sample_4_S6_indels_filtered_lifted_annotated_corrected.vep.vcf"
    "Exome-seq_Sample_5_S7_snvs_filtered_lifted_annotated_corrected.vep.vcf"
    "Exome-seq_Sample_5_S7_indels_filtered_lifted_annotated_corrected.vep.vcf"
    "Exome-seq_Sample_6_S8_snvs_filtered_lifted_annotated_corrected.vep.vcf"
    "Exome-seq_Sample_6_S8_indels_filtered_lifted_annotated_corrected.vep.vcf"
    "Exome-seq_Sample_7_S9_snvs_filtered_lifted_annotated_corrected.vep.vcf"
    "Exome-seq_Sample_7_S9_indels_filtered_lifted_annotated_corrected.vep.vcf"
    "Exome-seq_Sample_8_S10_snvs_filtered_lifted_annotated_corrected.vep.vcf"
    "Exome-seq_Sample_8_S10_indels_filtered_lifted_annotated_corrected.vep.vcf"
    "Exome-seq_Sample_9_S11_snvs_filtered_lifted_annotated_corrected.vep.vcf"
    "Exome-seq_Sample_9_S11_indels_filtered_lifted_annotated_corrected.vep.vcf"
    "Exome-seq_Sample_10_S12_snvs_filtered_lifted_annotated_corrected.vep.vcf"
    "Exome-seq_Sample_10_S12_indels_filtered_lifted_annotated_corrected.vep.vcf"
    "Exome-seq_Sample_11_S13_snvs_filtered_lifted_annotated_corrected.vep.vcf"
    "Exome-seq_Sample_11_S13_indels_filtered_lifted_annotated_corrected.vep.vcf"
    "Exome-seq_Sample_12_S14_snvs_filtered_lifted_annotated_corrected.vep.vcf"
    "Exome-seq_Sample_12_S14_indels_filtered_lifted_annotated_corrected.vep.vcf"
)

# Loop through each VCF file and annotate it with the dummy variable
for vcf_file in "${vcf_files[@]}"; do
    input_vcf="${BASE_DIR}/${vcf_file}"
    output_vcf="${BASE_DIR}/$(basename ${vcf_file} .vep.vcf)_with_dummy.vep.vcf"
    vcf-genotype-annotator "$input_vcf" DUMMY 0/1 -o "$output_vcf"
    echo "Annotated $input_vcf and saved to $output_vcf"
done

echo "All VCF files have been processed."
