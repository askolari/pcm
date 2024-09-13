import csv
import os

input_directory = '/mnt/c/Users/agsko/dev/pcm/DeepNeo'
output_directory = '/mnt/c/Users/agsko/dev/pcm/DeepNeo'
snv_prefix = 'S'
indel_prefix = 'S'
snv_suffix = '_snvs.combined.vep.vcf'
indel_suffix = '_indels.combined.vep.vcf'

header = [
    "Chromosome", "Position", "ID", "Ref", "Alt", "Quality", "Filter",
    "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type",
    "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position",
    "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation",
    "DISTANCE", "STRAND", "FLAGS", "SYMBOL_SOURCE", "HGNC_ID", "WildtypeProtein",
    "DownstreamProtein", "ProteinLengthChange"
]

def convert_vcf_to_csv(input_vcf, output_csv):
    with open(input_vcf, 'r') as infile, open(output_csv, 'w', newline='') as outfile:
        reader = (line for line in infile if not line.startswith('##'))
        reader = csv.reader(reader, delimiter='\t')

        writer = csv.writer(outfile)
        writer.writerow(header)

        for row in reader:
            if row[0].startswith('#'):
                continue  # Skip the header line

            chrom, pos, id_, ref, alt, qual, filter_, info, format_, normal, tumor = row
            info_fields = dict(field.split('=') for field in info.split(';') if '=' in field)

            if 'CSQ' not in info_fields:
                continue  # Skip lines without CSQ field

            csq_values = info_fields['CSQ'].split(',')
            for csq in csq_values:
                csq_data = csq.split('|')
                output_row = [chrom, pos, id_, ref, alt, qual, filter_] + csq_data
                writer.writerow(output_row)

    print(f"Converted {input_vcf} to {output_csv}")

for i in range(1, 16):
    snv_input_vcf = os.path.join(input_directory, f"{snv_prefix}{i}{snv_suffix}")
    snv_output_csv = os.path.join(output_directory, f"{snv_prefix}{i}{snv_suffix.replace('.vcf', '.csv')}")

    indel_input_vcf = os.path.join(input_directory, f"{indel_prefix}{i}{indel_suffix}")
    indel_output_csv = os.path.join(output_directory, f"{indel_prefix}{i}{indel_suffix.replace('.vcf', '.csv')}")

    if os.path.exists(snv_input_vcf):
        convert_vcf_to_csv(snv_input_vcf, snv_output_csv)
    else:
        print(f"File not found: {snv_input_vcf}")

    if os.path.exists(indel_input_vcf):
        convert_vcf_to_csv(indel_input_vcf, indel_output_csv)
    else:
        print(f"File not found: {indel_input_vcf}")

print("All conversions complete.")
