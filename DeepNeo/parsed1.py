import csv

input_vcf = '/mnt/c/Users/agsko/dev/pcm/DeepNeo/S1_snvs.combined.vep.vcf'
output_csv = '/mnt/c/Users/agsko/dev/pcm/DeepNeo/S1_snvs.combined.vep.parsed.csv'

header = [
    "Chromosome", "Position", "ID", "Ref", "Alt", "Quality", "Filter",
    "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type",
    "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position",
    "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation",
    "DISTANCE", "STRAND", "FLAGS", "SYMBOL_SOURCE", "HGNC_ID", "WildtypeProtein",
    "DownstreamProtein", "ProteinLengthChange"
]

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

print("VCF to CSV conversion complete.")
