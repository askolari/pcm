#!/bin/bash

# List of VCF files to process
vcf_files=(
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_1_S3_snvs_filtered_lifted_annotated_corrected_with_dummy.vep.vcf"
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_1_S3_indels_filtered_lifted_annotated_corrected_with_dummy.vep.vcf"
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_2_S4_snvs_filtered_lifted_annotated_corrected_with_dummy.vep.vcf"
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_2_S4_indels_filtered_lifted_annotated_corrected_with_dummy.vep.vcf"
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_3_S5_snvs_filtered_lifted_annotated_corrected_with_dummy.vep.vcf"
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_3_S5_indels_filtered_lifted_annotated_corrected_with_dummy.vep.vcf"
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_4_S6_snvs_filtered_lifted_annotated_corrected_with_dummy.vep.vcf"
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_4_S6_indels_filtered_lifted_annotated_corrected_with_dummy.vep.vcf"
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_5_S7_snvs_filtered_lifted_annotated_corrected_with_dummy.vep.vcf"
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_5_S7_indels_filtered_lifted_annotated_corrected_with_dummy.vep.vcf"
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_6_S8_snvs_filtered_lifted_annotated_corrected_with_dummy.vep.vcf"
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_6_S8_indels_filtered_lifted_annotated_corrected_with_dummy.vep.vcf"
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_7_S9_snvs_filtered_lifted_annotated_corrected_with_dummy.vep.vcf"
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_7_S9_indels_filtered_lifted_annotated_corrected_with_dummy.vep.vcf"
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_8_S10_snvs_filtered_lifted_annotated_corrected_with_dummy.vep.vcf"
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_8_S10_indels_filtered_lifted_annotated_corrected_with_dummy.vep.vcf"
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_9_S11_snvs_filtered_lifted_annotated_corrected_with_dummy.vep.vcf"
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_9_S11_indels_filtered_lifted_annotated_corrected_with_dummy.vep.vcf"
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_10_S12_snvs_filtered_lifted_annotated_corrected_with_dummy.vep.vcf"
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_10_S12_indels_filtered_lifted_annotated_corrected_with_dummy.vep.vcf"
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_11_S13_snvs_filtered_lifted_annotated_corrected_with_dummy.vep.vcf"
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_11_S13_indels_filtered_lifted_annotated_corrected_with_dummy.vep.vcf"
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_12_S14_snvs_filtered_lifted_annotated_corrected_with_dummy.vep.vcf"
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_12_S14_indels_filtered_lifted_annotated_corrected_with_dummy.vep.vcf"
)

# Generate mutant peptides for each VCF file
for vcf_file in "${vcf_files[@]}"
do
  output_file_class1="${vcf_file%.vep.vcf}_class1_mutant_peptides.fasta"
  output_file_class2="${vcf_file%.vep.vcf}_class2_mutant_peptides.fasta"

  # Generate Class I peptides
  pvacseq generate_protein_fasta $vcf_file 10 $output_file_class1 --mutant-only --sample-name DUMMY

  # Generate Class II peptides
  pvacseq generate_protein_fasta $vcf_file 15 $output_file_class2 --mutant-only --sample-name DUMMY
done

pvacseq generate_protein_fasta /mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_13_S15_indels_filtered_lifted_annotated_corrected_with_dummy.vep.vcf 10 /mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_13_S15_indels_filtered_lifted_annotated_corrected_with_dummy_class1_with_mutant_peptides.fasta --mutant-only --sample-name DUMMY
pvacseq generate_protein_fasta /mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_13_S15_indels_filtered_lifted_annotated_corrected_with_dummy.vep.vcf 15 /mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_13_S15_indels_filtered_lifted_annotated_corrected_with_dummy_class2_with_mutant_peptides.fasta --mutant-only --sample-name DUMMY
