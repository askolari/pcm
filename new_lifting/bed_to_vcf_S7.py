
import pandas as pd
import os

# Define the paths for the specific files
bed_file = "/mnt/c/Users/agsko/dev/pcm/new_lifting/Exome-seq_Sample_2_S4_indels_filtered_simplified_vcf_lifted.bed"
original_vcf_file = "/mnt/c/Users/agsko/dev/pcm/new_lifting/Exome-seq_Sample_2_S4_indels_filtered.vcf"
output_vcf_file = "/mnt/c/Users/agsko/dev/pcm/new_lifting/Exome-seq_Sample_2_S4_indels_filtered_lifted.vcf"

try:
    # Load the simplified lifted BED file
    lifted_bed = pd.read_csv(bed_file, sep='\t', header=None)
    lifted_bed.columns = ['CHROM', 'POS', 'END', 'ID']

    # Adjust positions: BED files are 0-based while VCF files are 1-based, so adjust the positions accordingly
    lifted_bed['POS'] = lifted_bed['POS'] + 1

    # Load the original VCF data (excluding the header)
    original_vcf = pd.read_csv(original_vcf_file, sep='\t', comment='#', header=None)
    original_vcf.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR']

    # Ensure both dataframes have the same length
    if len(original_vcf) != len(lifted_bed):
        raise ValueError("The original VCF and lifted BED files do not have the same number of entries.")

    # Merge the coordinates from the lifted BED into the original VCF
    merged_vcf = original_vcf.copy()
    merged_vcf['CHROM'] = lifted_bed['CHROM']
    merged_vcf['POS'] = lifted_bed['POS']

    # Load the original VCF header
    with open(original_vcf_file, 'r') as file:
        vcf_header = []
        for line in file:
            if line.startswith("#"):
                vcf_header.append(line.strip())
            else:
                break

    # Save the VCF file with the updated coordinates
    with open(output_vcf_file, 'w') as file:
        # Write the VCF header
        for line in vcf_header:
            file.write(line + '\n')
        # Write the updated VCF data
        merged_vcf.to_csv(file, sep='\t', index=False, header=False)

    print(f"Lifted VCF file created: {output_vcf_file}")
except Exception as e:
    print(f"An error occurred: {e}")




perl /mnt/c/Users/agsko/dev/pcm/ensembl-vep/vep \
  -i /mnt/c/Users/agsko/dev/pcm/new_lifting/Exome-seq_Sample_13_S15_indels_filtered_lifted.vcf \
  -o /mnt/c/Users/agsko/dev/pcm/new_lifting/Exome-seq_Sample_13_S15_indel_filtered_lifted_annotated.vep.vcf \
  --vcf --symbol --terms SO --tsl --biotype --hgvs \
  --fasta /mnt/c/Users/agsko/dev/pcm/ensembl-vep/fasta/Mus_musculus.GRCm38.dna.primary_assembly.fa \
  --offline --cache \
  --dir_cache /mnt/c/Users/agsko/dev/pcm/ensembl-vep/cache \
  --plugin Frameshift --plugin Wildtype --plugin ProteinSeqs --plugin LoFtool --plugin NMD \
  --dir_plugins /mnt/c/Users/agsko/dev/pcm/ensembl-vep/Plugins \
  --species mus_musculus \
  --pick --transcript_version