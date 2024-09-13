import pandas as pd
import os

# Function to parse the INFO column manually
def parse_info_column(info_str):
    info_dict = {}
    for entry in info_str.split(';'):
        if '=' in entry:
            key, value = entry.split('=', 1)
            if key == 'CSQ':
                info_dict[key] = value.split(',')
            else:
                info_dict[key] = value
        else:
            info_dict[entry] = True
    return info_dict

# Function to parse individual CSQ entries
def parse_csq(csq_str):
    csq_fields = [
        "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", 
        "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", 
        "Protein_position", "Amino_acids", "Codons", "Existing_variation", "DISTANCE", 
        "STRAND", "FLAGS", "VARIANT_CLASS", "SYMBOL_SOURCE", "HGNC_ID", "CANONICAL", 
        "MANE_SELECT", "MANE_PLUS_CLINICAL", "TSL", "APPRIS", "CCDS", "ENSP", "SWISSPROT", 
        "TREMBL", "UNIPARC", "UNIPROT_ISOFORM", "GENE_PHENO", "SIFT", "DOMAINS", "miRNA", 
        "HGVS_OFFSET", "AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF", "gnomADe_AF", 
        "gnomADe_AFR_AF", "gnomADe_AMR_AF", "gnomADe_ASJ_AF", "gnomADe_EAS_AF", "gnomADe_FIN_AF", 
        "gnomADe_NFE_AF", "gnomADe_OTH_AF", "gnomADe_SAS_AF", "gnomADg_AF", "gnomADg_AFR_AF", 
        "gnomADg_AMI_AF", "gnomADg_AMR_AF", "gnomADg_ASJ_AF", "gnomADg_EAS_AF", "gnomADg_FIN_AF", 
        "gnomADg_MID_AF", "gnomADg_NFE_AF", "gnomADg_OTH_AF", "gnomADg_SAS_AF", "MAX_AF", 
        "MAX_AF_POPS", "CLIN_SIG", "SOMATIC", "PHENO", "PUBMED", "MOTIF_NAME", "MOTIF_POS", 
        "HIGH_INF_POS", "MOTIF_SCORE_CHANGE", "TRANSCRIPTION_FACTORS"
    ]
    csq_values = csq_str.split('|')
    return dict(zip(csq_fields, csq_values))

# Function to expand the CSQ field into individual columns
def expand_csq(df):
    csq_expanded = df['CSQ'].explode().apply(pd.Series)
    return df.drop(columns=['CSQ']).join(csq_expanded, rsuffix='_CSQ')

# Function to process each VCF file
def process_file(file_path):
    try:
        df = pd.read_csv(file_path, sep='\t', comment='#', header=None)
    except FileNotFoundError:
        print(f"File {file_path} not found, skipping...")
        return
    except pd.errors.EmptyDataError:
        print(f"No columns to parse from {file_path}, skipping...")
        return

    # Assuming INFO is the 8th column in the VCF
    df['INFO_parsed'] = df.iloc[:, 7].apply(parse_info_column)
    
    # Expand the INFO column into separate columns
    info_df = df['INFO_parsed'].apply(pd.Series)
    
    # Further expand the CSQ field into individual columns
    if 'CSQ' in info_df.columns:
        csq_df = pd.DataFrame([parse_csq(csq) for csq in info_df['CSQ'].explode()])
        info_df = info_df.drop(columns=['CSQ']).join(csq_df, rsuffix='_CSQ')
    
    # Concatenate the original DataFrame with the expanded INFO DataFrame
    combined_df = pd.concat([df.drop(columns=['INFO_parsed']), info_df], axis=1)
    
    # Save the processed DataFrame to a new CSV file
    new_file_path = file_path.replace('.vcf', '_parsed.csv')
    combined_df.to_csv(new_file_path, index=False)
    print(f"Processed file saved as: {new_file_path}")

# List of VCF file paths
vcf_files = [
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_12_S14_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_12_S14_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_13_S15_indel_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_13_S15_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Control_1_S1_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Control_2_S2_snvs_filtered_lifted_annotated.vep.vcf",
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
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_9_S11_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_10_S12_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_10_S12_snvs_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_11_S13_indels_filtered_lifted_annotated.vep.vcf",
    "/mnt/c/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_11_S13_snvs_filtered_lifted_annotated.vep.vcf",
]

# Process each VCF file to create parsed CSV files
for vcf_file in vcf_files:
    process_file(vcf_file)
