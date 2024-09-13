import pandas as pd
import os

# Load the CSV files
combined_final_peptides_path = '/mnt/c/Users/agsko/dev/pcm/predictions/combined_final_peptides.csv'
all_high_affinity_peptides_path = '/mnt/c/Users/agsko/dev/pcm/predictions/all_high_affinity_peptides_combined.csv'

combined_final_peptides = pd.read_csv(combined_final_peptides_path)
all_high_affinity_peptides = pd.read_csv(all_high_affinity_peptides_path)

# Extract the relevant columns from the combined_final_peptides dataframe
relevant_columns = ['Gene', 'Transcript', 'Transcript_Version', 'Type', 'Position']

# Create new columns in the all_high_affinity_peptides DataFrame
all_high_affinity_peptides['Matching_Sequence_1'] = ''
all_high_affinity_peptides['Gene_1'] = ''
all_high_affinity_peptides['Matching_Sequence_2'] = ''
all_high_affinity_peptides['Gene_2'] = ''
# Add more columns if needed

# Loop over each row in the all_high_affinity_peptides DataFrame
for index, row in all_high_affinity_peptides.iterrows():
    peptide = row['peptide']
    sample_id = row['Sample_ID']
    
    # Find matching sequences in the combined_final_peptides DataFrame
    matches = combined_final_peptides[(combined_final_peptides['Sequence'].str.contains(peptide)) & 
                                      (combined_final_peptides['Sample_ID'] == sample_id)]
    
    # Populate the new columns with matching sequences and corresponding genes
    if not matches.empty:
        match_count = 1
        for _, match_row in matches.iterrows():
            all_high_affinity_peptides.at[index, f'Matching_Sequence_{match_count}'] = match_row['Sequence']
            
            # Add the relevant columns
            for col in relevant_columns:
                all_high_affinity_peptides.at[index, f'{col}_{match_count}'] = match_row[col]
            
            match_count += 1


# Save the updated DataFrame to a new CSV file
output_path = '/mnt/c/Users/agsko/dev/pcm/predictions/all_high_affinity_pred_class1.csv'
all_high_affinity_peptides.to_csv(output_path, index=False)

print(f"Updated data saved to {output_path}")

import pandas as pd
import os

# Load the CSV files
combined_final_peptides_path = '/mnt/c/Users/agsko/dev/pcm/predictions/combined_final_peptides.csv'
all_high_affinity_peptides_path = '/mnt/c/Users/agsko/dev/pcm/predictions/all_high_affinity_class2_peptides_combined.csv'

combined_final_peptides = pd.read_csv(combined_final_peptides_path)
all_high_affinity_peptides = pd.read_csv(all_high_affinity_peptides_path)

# Extract the relevant columns from the combined_final_peptides dataframe
relevant_columns = ['Gene', 'Transcript', 'Transcript_Version', 'Type', 'Position']

# Create new columns in the all_high_affinity_peptides DataFrame
all_high_affinity_peptides['Matching_Sequence_1'] = ''
all_high_affinity_peptides['Gene_1'] = ''
all_high_affinity_peptides['Matching_Sequence_2'] = ''
all_high_affinity_peptides['Gene_2'] = ''
# Add more columns if needed

# Loop over each row in the all_high_affinity_peptides DataFrame
for index, row in all_high_affinity_peptides.iterrows():
    peptide = row['peptide']
    sample_id = row['Sample_ID']
    
    # Find matching sequences in the combined_final_peptides DataFrame
    matches = combined_final_peptides[(combined_final_peptides['Sequence'].str.contains(peptide)) & 
                                      (combined_final_peptides['Sample_ID'] == sample_id)]
    
    # Populate the new columns with matching sequences and corresponding genes
    if not matches.empty:
        match_count = 1
        for _, match_row in matches.iterrows():
            all_high_affinity_peptides.at[index, f'Matching_Sequence_{match_count}'] = match_row['Sequence']
            
            # Add the relevant columns
            for col in relevant_columns:
                all_high_affinity_peptides.at[index, f'{col}_{match_count}'] = match_row[col]
            
            match_count += 1


# Save the updated DataFrame to a new CSV file
output_path = '/mnt/c/Users/agsko/dev/pcm/predictions/all_high_affinity_pred_class2.csv'
all_high_affinity_peptides.to_csv(output_path, index=False)

print(f"Updated data saved to {output_path}")
