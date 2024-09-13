import os

# Define the input and output file paths
input_vcf = "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_13_S15_snvs_filtered_lifted_annotated.vep.vcf"
output_vcf = input_vcf.replace(".vep.vcf", "_corrected.vep.vcf")

# Correct the QUAL field in the VCF
with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
    for line in infile:
        if line.startswith('#'):
            outfile.write(line)
        else:
            fields = line.strip().split('\t')
            if not fields[5].replace('.', '', 1).isdigit():
                fields[5] = '.'
            outfile.write('\t'.join(fields) + '\n')

print(f"Corrected VCF file created: {output_vcf}")
