import os

# List of your VCF file paths
vcf_files = [
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_10_S12_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_10_S12_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_11_S13_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_11_S13_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_12_S14_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_12_S14_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_1_S3_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_1_S3_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_2_S4_indel_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_2_S4_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_3_S5_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_3_S5_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_4_S6_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_4_S6_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_5_S7_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_5_S7_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_6_S8_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_6_S8_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_7_S9_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_7_S9_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_8_S10_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_8_S10_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_9_S11_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_9_S11_snvs_filtered_lifted_annotated.vep.vcf"
]

# Loop through each file, correct it, and save it
for vcf_file in vcf_files:
    output_vcf_file = vcf_file.replace(".vep.vcf", "_corrected.vep.vcf")
    
    with open(vcf_file, 'r') as infile, open(output_vcf_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
            else:
                fields = line.strip().split('\t')
                if not fields[5].replace('.', '', 1).isdigit():
                    fields[5] = '.'
                outfile.write('\t'.join(fields) + '\n')

    print(f"VCF file {vcf_file} has been corrected and saved as {output_vcf_file}.")
