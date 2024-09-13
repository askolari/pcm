import os

def vcf_to_bed_simplified(vcf_file, bed_file):
    with open(vcf_file, 'r') as vcf, open(bed_file, 'w') as bed:
        for line in vcf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            start = pos - 1  # BED format uses 0-based start position
            ref = fields[3]  # Assuming this to be the ID or key identifier for simplification
            bed.write(f"{chrom}\t{start}\t{pos}\t{ref}\n")

# Define directory and file paths
vcf_dir = '/mnt/c/Users/agsko/dev/pcm/new_lifting/'
vcf_files = [
    "Exome-seq_Control_1_S1_indels_filtered.vcf",
    "Exome-seq_Control_1_S1_snvs_filtered.vcf",
    "Exome-seq_Control_2_S2_indels_filtered.vcf",
    "Exome-seq_Control_2_S2_snvs_filtered.vcf",
    "Exome-seq_Sample_1_S3_indels_filtered.vcf",
    "Exome-seq_Sample_1_S3_snvs_filtered.vcf",
    "Exome-seq_Sample_2_S4_indels_filtered.vcf",
    "Exome-seq_Sample_2_S4_snvs_filtered.vcf",
    "Exome-seq_Sample_3_S5_indels_filtered.vcf",
    "Exome-seq_Sample_3_S5_snvs_filtered.vcf",
    "Exome-seq_Sample_4_S6_indels_filtered.vcf",
    "Exome-seq_Sample_4_S6_snvs_filtered.vcf",
    "Exome-seq_Sample_5_S7_indels_filtered.vcf",
    "Exome-seq_Sample_5_S7_snvs_filtered.vcf",
    "Exome-seq_Sample_6_S8_indels_filtered.vcf",
    "Exome-seq_Sample_6_S8_snvs_filtered.vcf",
    "Exome-seq_Sample_7_S9_indels_filtered.vcf",
    "Exome-seq_Sample_7_S9_snvs_filtered.vcf",
    "Exome-seq_Sample_8_S10_indels_filtered.vcf",
    "Exome-seq_Sample_8_S10_snvs_filtered.vcf",
    "Exome-seq_Sample_9_S11_indels_filtered.vcf",
    "Exome-seq_Sample_9_S11_snvs_filtered.vcf",
    "Exome-seq_Sample_10_S12_indels_filtered.vcf",
    "Exome-seq_Sample_10_S12_snvs_filtered.vcf",
    "Exome-seq_Sample_11_S13_indels_filtered.vcf",
    "Exome-seq_Sample_11_S13_snvs_filtered.vcf",
    "Exome-seq_Sample_12_S14_indels_filtered.vcf",
    "Exome-seq_Sample_12_S14_snvs_filtered.vcf",
    "Exome-seq_Sample_13_S15_indels_filtered.vcf",
    "Exome-seq_Sample_13_S15_snvs_filtered.vcf"
]

# Process each VCF file
for vcf_file in vcf_files:
    vcf_path = os.path.join(vcf_dir, vcf_file)
    bed_file = vcf_path.replace('.vcf', '_simplified.bed')
    vcf_to_bed_simplified(vcf_path, bed_file)
    print(f"Converted {vcf_file} to BED format.")
