def vcf_to_bed_extended(vcf_file, bed_file):
    with open(vcf_file, 'r') as vcf, open(bed_file, 'w') as bed:
        for line in vcf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            qual = fields[5]
            filt = fields[6]
            info = fields[7]
            format_field = fields[8]
            normal_sample = fields[9]
            tumor_sample = fields[10]
            bed.write(f"{chrom}\t{int(pos)-1}\t{pos}\t{ref}\t{alt}\t{qual}\t{filt}\t{info}\t{format_field}\t{normal_sample}\t{tumor_sample}\n")

# Usage
vcf_file = '/mnt/c/Users/agsko/dev/pcm/Exome-seq_Control_2_S2_snvs_filtered.vcf'
bed_file = '/mnt/c/Users/agsko/dev/pcm/Exome-seq_Control_2_S2_snvs_filtered_extended.bed'
vcf_to_bed_extended(vcf_file, bed_file)
