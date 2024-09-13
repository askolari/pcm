import pandas as pd

# Paths to your files
vcf_file = '/mnt/c/Users/agsko/dev/pcm/Exome-seq_Sample_13_S15_snvs_filtered_lifted_annotated.vep.vcf'
output_vcf_file = '/mnt/c/Users/agsko/dev/pcm/Exome-seq_Sample_13_S15_snvs_filtered_lifted_annotated_corrected.vep.vcf'

# Open and process the VCF file
with open(vcf_file, 'r') as infile, open(output_vcf_file, 'w') as outfile:
    for line in infile:
        if line.startswith('#'):
            outfile.write(line)
        else:
            fields = line.strip().split('\t')
            if not fields[5].replace('.', '', 1).isdigit():
                fields[5] = '.'
            outfile.write('\t'.join(fields) + '\n')

print("VCF file has been corrected and saved.")
