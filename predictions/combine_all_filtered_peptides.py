import pandas as pd

# Define the file paths for the datasets
file_paths = [
    "/mnt/c/Users/agsko/dev/pcm/predictions/X101_matched_class1.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X88_matched_class1.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X89_matched_class1.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X90_matched_class1.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X92_matched_class1.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X93_matched_class1.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X94_matched_class1.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X97_matched_class1.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X99_matched_class1.csv"
]

# Initialize an empty list to store dataframes
dataframes = []

# Loop through the file paths, read each CSV file, and append it to the list of dataframes
for file_path in file_paths:
    df = pd.read_csv(file_path)
    dataframes.append(df)

# Concatenate all dataframes into a single dataframe
combined_df = pd.concat(dataframes, ignore_index=True)

# Save the combined dataframe to a new CSV file
output_file = "/mnt/c/Users/agsko/dev/pcm/predictions/all_peptides_class1_filtered_by_expression.csv"
combined_df.to_csv(output_file, index=False)

print(f"Combined data saved to {output_file}")

import pandas as pd

# Define the file paths for class 2 datasets
file_paths_class2 = [
    "/mnt/c/Users/agsko/dev/pcm/predictions/X101_matched_class2.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X88_matched_class2.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X89_matched_class2.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X90_matched_class2.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X92_matched_class2.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X93_matched_class2.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X94_matched_class2.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X97_matched_class2.csv",
    "/mnt/c/Users/agsko/dev/pcm/predictions/X99_matched_class2.csv"
]

# Initialize an empty list to store dataframes
dataframes_class2 = []

# Loop through the file paths, read each CSV file, and append it to the list of dataframes
for file_path in file_paths_class2:
    df_class2 = pd.read_csv(file_path)
    dataframes_class2.append(df_class2)

# Concatenate all dataframes into a single dataframe for class 2
combined_df_class2 = pd.concat(dataframes_class2, ignore_index=True)

# Save the combined dataframe to a new CSV file for class 2
output_file_class2 = "/mnt/c/Users/agsko/dev/pcm/predictions/all_peptides_class2_filtered_by_expression.csv"
combined_df_class2.to_csv(output_file_class2, index=False)

print(f"Combined data for class 2 saved to {output_file_class2}")

import pandas as pd

# Load the class 1 and class 2 datasets
class1_file = "/mnt/c/Users/agsko/dev/pcm/predictions/all_peptides_class1_filtered_by_expression.csv"
class2_file = "/mnt/c/Users/agsko/dev/pcm/predictions/all_peptides_class2_filtered_by_expression.csv"

# Load the data
class1_df = pd.read_csv(class1_file)
class2_df = pd.read_csv(class2_file)

# Remove duplicates based on 'peptide', 'netmhcpan_ba.IC50', and 'Sample.ID' columns
class1_unique_df = class1_df.drop_duplicates(subset=['peptide', 'netmhcpan_ba.IC50', 'Sample.ID'])
class2_unique_df = class2_df.drop_duplicates(subset=['peptide', 'ic50', 'Sample.ID'])

# Save the new files as *_final.csv
class1_output_file = "/mnt/c/Users/agsko/dev/pcm/predictions/all_peptides_class1_filtered_by_expression_final.csv"
class2_output_file = "/mnt/c/Users/agsko/dev/pcm/predictions/all_peptides_class2_filtered_by_expression_final.csv"

class1_unique_df.to_csv(class1_output_file, index=False)
class2_unique_df.to_csv(class2_output_file, index=False)

print(f"Final class 1 data saved to {class1_output_file}")
print(f"Final class 2 data saved to {class2_output_file}")

import pandas as pd

# Load the final class 1 and class 2 datasets
class1_file = "/mnt/c/Users/agsko/dev/pcm/predictions/all_peptides_class1_filtered_by_expression_final.csv"
class2_file = "/mnt/c/Users/agsko/dev/pcm/predictions/all_peptides_class2_filtered_by_expression_final.csv"

# Load the data
class1_df = pd.read_csv(class1_file)
class2_df = pd.read_csv(class2_file)

# Check for duplicated peptides in Class 1
duplicated_peptides_class1 = class1_df.groupby('peptide').filter(lambda x: x['Sample.ID'].nunique() > 1)

# Check for duplicated peptides in Class 2
duplicated_peptides_class2 = class2_df.groupby('peptide').filter(lambda x: x['Sample.ID'].nunique() > 1)

# Print and save the results for Class 1
if not duplicated_peptides_class1.empty:
    print("Peptides found in more than one sample (Class 1):")
    print(duplicated_peptides_class1)
    
    # Save duplicated peptides for Class 1
    output_class1_file = "/mnt/c/Users/agsko/dev/pcm/predictions/peptides_in_multiple_samples_class1_final.csv"
    duplicated_peptides_class1.to_csv(output_class1_file, index=False)
    print(f"Duplicated peptides for Class 1 saved to {output_class1_file}")
else:
    print("No peptides found in more than one sample for Class 1.")

# Print and save the results for Class 2
if not duplicated_peptides_class2.empty:
    print("Peptides found in more than one sample (Class 2):")
    print(duplicated_peptides_class2)
    
    # Save duplicated peptides for Class 2
    output_class2_file = "/mnt/c/Users/agsko/dev/pcm/predictions/peptides_in_multiple_samples_class2_final.csv"
    duplicated_peptides_class2.to_csv(output_class2_file, index=False)
    print(f"Duplicated peptides for Class 2 saved to {output_class2_file}")
else:
    print("No peptides found in more than one sample for Class 2.")
