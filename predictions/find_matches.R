
# Load necessary libraries
library(dplyr)
library(tidyr)

# Function to parse the 'Header' column and split into new columns
parse_fasta_data <- function(csv_file_path, output_csv_path) {
  # Read the CSV file
  df <- read.csv(csv_file_path, stringsAsFactors = FALSE)
  
  # Split the 'Header' column into multiple columns based on the period (.)
  df_parsed <- df %>%
    separate(Header, into = c("Prefix", "ID1", "Gene", "Transcript", "Version", "Type", "Mutation"), sep = "\\.") 
  
  # Save the parsed data frame to a new CSV file
  write.csv(df_parsed, output_csv_path, row.names = FALSE)
  cat("Parsed and saved:", output_csv_path, "\n")
}


file_paths <- list(
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_1_S3_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_1_S3_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_1_S3_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_1_S3_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_2_S4_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_2_S4_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_2_S4_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_2_S4_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_3_S5_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_3_S5_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_3_S5_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_3_S5_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_4_S6_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_4_S6_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_4_S6_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_4_S6_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_5_S7_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_5_S7_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_5_S7_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_5_S7_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_6_S8_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_6_S8_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_6_S8_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_6_S8_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_7_S9_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_7_S9_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_7_S9_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_7_S9_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_8_S10_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_8_S10_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_8_S10_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_8_S10_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_9_S11_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_9_S11_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_9_S11_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_9_S11_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_10_S12_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_10_S12_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_10_S12_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_10_S12_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_11_S13_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_11_S13_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_11_S13_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_11_S13_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_12_S14_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_12_S14_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_12_S14_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_12_S14_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_13_S15_indels_filtered_lifted_annotated_corrected_with_dummy_class1_with_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_13_S15_indels_filtered_lifted_annotated_corrected_with_dummy_class2_with_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_13_S15_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_with_mutant_peptides.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/Exome-seq_Sample_13_S15_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_with_mutant_peptides.csv"
)
# Directory to save the output CSV files
output_directory <- "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files"

# Create the output directory if it doesn't exist
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

# Loop through each file, parse the header, and save the result to a new file
for (file_path in file_paths) {
  # Generate the output file name
  csv_file_name <- sub(".csv$", "_parsed.csv", basename(file_path))
  output_csv_path <- file.path(output_directory, csv_file_name)
  
  # Parse and save the file
  parse_fasta_data(file_path, output_csv_path)
}


# List of input CSV file paths
file_paths <- list(
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_13_S15_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_with_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_13_S15_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_with_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_13_S15_indels_filtered_lifted_annotated_corrected_with_dummy_class2_with_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_13_S15_indels_filtered_lifted_annotated_corrected_with_dummy_class1_with_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_12_S14_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_12_S14_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_12_S14_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_12_S14_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_11_S13_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_11_S13_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_11_S13_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_11_S13_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_10_S12_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_10_S12_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_10_S12_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_10_S12_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_9_S11_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_9_S11_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_9_S11_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_9_S11_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_8_S10_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_8_S10_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_8_S10_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_8_S10_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_7_S9_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_7_S9_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_7_S9_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_7_S9_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_6_S8_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_6_S8_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_6_S8_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_6_S8_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_5_S7_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_5_S7_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_5_S7_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_5_S7_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_4_S6_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_4_S6_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_4_S6_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_4_S6_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_3_S5_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_3_S5_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_3_S5_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_3_S5_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_2_S4_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_2_S4_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_2_S4_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_2_S4_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_1_S3_snvs_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_1_S3_snvs_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_1_S3_indels_filtered_lifted_annotated_corrected_with_dummy_class2_mutant_peptides_parsed.csv",
  "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/Exome-seq_Sample_1_S3_indels_filtered_lifted_annotated_corrected_with_dummy_class1_mutant_peptides_parsed.csv"
)

# Initialize an empty list to store each CSV data
all_data <- list()

# Loop through each file path, read the CSV, and add the sample information
for (file_path in file_paths) {
  
  # Extract the full file name (basename of the path)
  file_name <- basename(file_path)
  
  # Try to extract the sample ID (e.g., S3). If it doesn't work, keep the full file name.
  # Use regex to extract the part like "S3", assuming it is always after "_Sample_" and before "_"
  sample_id <- sub(".*_Sample_([A-Za-z0-9]+)_.*", "\\1", file_name)
  
  # If sample_id is still the same as the file name, it means extraction failed. Use the full file name instead.
  if (sample_id == file_name) {
    sample_id <- file_name
  }
  
  # Read the CSV file
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Add columns for the extracted sample ID and mutation type (snv or indel)
  mutation_type <- ifelse(grepl("snvs", file_name), "snv", "indel")
  data$Sample <- sample_id
  data$Mutation_Type <- mutation_type
  
  # Append the data to the list
  all_data[[file_name]] <- data
}

# Combine all the data into one dataframe
combined_data <- bind_rows(all_data)

# Save the combined data into a new CSV file
write.csv(combined_data, "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/combined_data.csv", row.names = FALSE)

cat("Combined CSV has been saved successfully!\n")

# List of files and corresponding sample IDs
file_paths <- c(
  "C:/Users/agsko/dev/pcm/predictions/S15_class1_predictions.csv",
  "C:/Users/agsko/dev/pcm/predictions/S3_class1_predictions.csv",
  "C:/Users/agsko/dev/pcm/predictions/S4_class1_predictions.csv",
  "C:/Users/agsko/dev/pcm/predictions/S5_class1_predictions.csv",
  "C:/Users/agsko/dev/pcm/predictions/S6_class1_predictions.csv",
  "C:/Users/agsko/dev/pcm/predictions/S7_class1_predictions.csv",
  "C:/Users/agsko/dev/pcm/predictions/S8_class1_predictions.csv",
  "C:/Users/agsko/dev/pcm/predictions/S9_class1_predictions.csv",
  "C:/Users/agsko/dev/pcm/predictions/S10_class1_predictions.csv",
  "C:/Users/agsko/dev/pcm/predictions/S11_class1_predictions.csv",
  "C:/Users/agsko/dev/pcm/predictions/S12_class1_predictions.csv",
  "C:/Users/agsko/dev/pcm/predictions/S13_class1_predictions.csv",
  "C:/Users/agsko/dev/pcm/predictions/S14_class1_predictions.csv"
)

sample_ids <- c("S15", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "S13", "S14")

# Initialize an empty list to store data frames
combined_data <- list()

# Loop through the files and add the Sample ID column
for (i in 1:length(file_paths)) {
  # Read the CSV file
  df <- read.csv(file_paths[i])
  
  # Add a new column for Sample ID
  df$Sample_ID <- sample_ids[i]
  
  # Append to the list
  combined_data[[i]] <- df
}

# Combine all data frames into one
final_combined_data <- do.call(rbind, combined_data)

# Save the combined data to a CSV file
write.csv(final_combined_data, "C:/Users/agsko/dev/pcm/predictions/combined_class1_predictions.csv", row.names = FALSE)

# Print the first few rows to check the output
head(final_combined_data)

# List of files and corresponding sample IDs for class2
file_paths <- c(
  "C:/Users/agsko/dev/pcm/predictions/S15_class2_predictions.csv",
  "C:/Users/agsko/dev/pcm/predictions/S3_class2_predictions.csv",
  "C:/Users/agsko/dev/pcm/predictions/S4_class2_predictions.csv",
  "C:/Users/agsko/dev/pcm/predictions/S5_class2_predictions.csv",
  "C:/Users/agsko/dev/pcm/predictions/S6_class2_predictions.csv",
  "C:/Users/agsko/dev/pcm/predictions/S7_class2_predictions.csv",
  "C:/Users/agsko/dev/pcm/predictions/S8_class2_predictions.csv",
  "C:/Users/agsko/dev/pcm/predictions/S9_class2_predictions.csv",
  "C:/Users/agsko/dev/pcm/predictions/S10_class2_predictions.csv",
  "C:/Users/agsko/dev/pcm/predictions/S11_class2_predictions.csv",
  "C:/Users/agsko/dev/pcm/predictions/S12_class2_predictions.csv",
  "C:/Users/agsko/dev/pcm/predictions/S13_class2_predictions.csv",
  "C:/Users/agsko/dev/pcm/predictions/S14_class2_predictions.csv"
)

sample_ids <- c("S15", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "S13", "S14")

# Initialize an empty list to store data frames
combined_data <- list()

# Loop through the files and add the Sample ID column
for (i in 1:length(file_paths)) {
  # Read the CSV file
  df <- read.csv(file_paths[i])
  
  # Add a new column for Sample ID
  df$Sample_ID <- sample_ids[i]
  
  # Append to the list
  combined_data[[i]] <- df
}

# Combine all data frames into one
final_combined_data <- do.call(rbind, combined_data)

# Save the combined data to a CSV file
write.csv(final_combined_data, "C:/Users/agsko/dev/pcm/predictions/combined_class2_predictions.csv", row.names = FALSE)

# Print the first few rows to check the output
head(final_combined_data)

# Load the combined class 1 and class 2 predictions
class1_predictions <- read.csv("C:/Users/agsko/dev/pcm/predictions/combined_class1_predictions.csv")
class2_predictions <- read.csv("C:/Users/agsko/dev/pcm/predictions/combined_class2_predictions.csv")

# Filter rows where netmhcpan_ba.IC50 < 500 for class 1
class1_filtered <- subset(class1_predictions, netmhcpan_ba.IC50 < 500)

# Filter rows where ic50 < 500 for class 2
class2_filtered <- subset(class2_predictions, ic50 < 500)

# Save the filtered results to new CSV files
write.csv(class1_filtered, "C:/Users/agsko/dev/pcm/predictions/combined_class1_predictions_IC500.csv", row.names = FALSE)
write.csv(class2_filtered, "C:/Users/agsko/dev/pcm/predictions/combined_class2_predictions_IC500.csv", row.names = FALSE)

# Print a few rows to check the output
head(class1_filtered)
head(class2_filtered)


# Define file paths
peptide_file_path <- "C:/Users/agsko/dev/pcm/predictions/combined_class1_predictions_IC500.csv"
combined_data_file_path <- "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/combined_data.csv"
output_file_path <- "C:/Users/agsko/dev/pcm/predictions/matched_class1_all_peptides.csv"

# Read the peptides file
peptides_data <- read.csv(peptide_file_path, stringsAsFactors = FALSE)

# Read the combined data with sequences
combined_data <- read.csv(combined_data_file_path, stringsAsFactors = FALSE)

# Initialize an empty data frame to store matched peptides
matched_data <- data.frame()

# Loop through each peptide and search for it as a substring in the combined_data sequences
for (i in 1:nrow(peptides_data)) {
  peptide <- peptides_data$peptide[i]
  
  # Check if the peptide exists as a substring in any sequence in combined_data
  matches <- combined_data[grep(peptide, combined_data$Sequence), ]
  
  # If matches are found, add the relevant information to the matched_data
  if (nrow(matches) > 0) {
    # Add the matched rows along with the corresponding peptide data
    for (j in 1:nrow(matches)) {
      matched_row <- cbind(peptides_data[i, ], matches[j, c("Gene","Transcript", "Type", "Mutation_Type", "Sample")])
      matched_data <- rbind(matched_data, matched_row)
    }
  }
}

# Save the matched data to a new CSV file
write.csv(matched_data, output_file_path, row.names = FALSE)

cat("Matched peptides file saved successfully!\n")

# Define the file path
file_path <- "C:/Users/agsko/dev/pcm/predictions/matched_class1_all_peptides.csv"

# Load the CSV file into a DataFrame
df <- read.csv(file_path)

# Check for duplicates based on all columns
duplicates <- df %>%
  filter(duplicated(.) | duplicated(., fromLast = TRUE))  # Find all duplicate rows

# Display duplicates
if (nrow(duplicates) > 0) {
  print("Found duplicates:")
  print(duplicates)
} else {
  print("No duplicates found.")
}

# Optionally, save the duplicates to a separate file
duplicates_file <- "C:/Users/agsko/dev/pcm/predictions/duplicate_class1_all_peptides.csv"
write.csv(duplicates, file = duplicates_file, row.names = FALSE)
print(paste("Duplicate rows saved to", duplicates_file))


# Define the file path
file_path <- "C:/Users/agsko/dev/pcm/predictions/matched_class1_all_peptides.csv"

# Load the CSV file into a DataFrame
df <- read.csv(file_path)

# Remove the second occurrence of duplicate rows (keep only the first occurrence)
df_clean <- df %>%
  distinct()  # This removes duplicate rows based on all columns

# Save the clean dataset to a new file
clean_file <- "C:/Users/agsko/dev/pcm/predictions/cleaned_class1_all_peptides.csv"
write.csv(df_clean, file = clean_file, row.names = FALSE)

# Load your CSV file
clean_file <- "C:/Users/agsko/dev/pcm/predictions/cleaned_class1_all_peptides.csv"
cleaned_peptides <- read.csv(clean_file)

# Add the Sample_Group column based on the Sample number
cleaned_peptides$Sample_Group <- with(cleaned_peptides, 
                                      ifelse(Sample %in% 1:4, "Ure_mock_IR", 
                                             ifelse(Sample %in% 5:8, "Ure_mock_IR_DMXAA", 
                                                    ifelse(Sample %in% 9:12, "Ure_IR", 
                                                           ifelse(Sample == 13, "PBS_IR_DMXAA", NA)))))

# Save the updated file under the same name
write.csv(cleaned_peptides, clean_file, row.names = FALSE)

# Print the first few rows to verify
head(cleaned_peptides)


# Print the number of rows before and after cleaning
print(paste("Original number of rows:", nrow(df)))
print(paste("Number of rows after removing duplicates:", nrow(df_clean)))

# Print a message
print(paste("Cleaned dataset saved to", clean_file))

# MANUALLY change Transcript to Feature column

# Load your CSV files
cleaned_peptides <- read.csv("C:/Users/agsko/dev/pcm/predictions/cleaned_class1_all_peptides.csv")
final_variants <- read.csv("C:/Users/agsko/dev/pcm/predictions/final_variants_mm10.csv")

# Remove rows from final_variants where Feature is empty
final_variants <- final_variants[final_variants$Feature != "", ]

# Function to perform the partial matching, safely handling special characters
find_partial_matches <- function(peptide_feature, variants_features) {
  # Escape any special characters in the peptide feature
  peptide_feature_escaped <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", peptide_feature)
  
  # Check if any variant feature contains the peptide feature as a substring
  matches <- grepl(peptide_feature_escaped, variants_features)
  return(matches)
}

# Initialize an empty list to store results
matched_results <- list()

# Loop through the cleaned_peptides features
for (i in 1:nrow(cleaned_peptides)) {
  peptide_feature <- cleaned_peptides$Feature[i]
  sample_group <- cleaned_peptides$Sample_Group[i]
  
  # Find matches in final_variants Feature column
  matches <- find_partial_matches(peptide_feature, final_variants$Feature)
  
  # If matches are found
  if (any(matches)) {
    # Subset final_variants based on the matches
    matched_variants <- final_variants[matches & final_variants$Sample_Group == sample_group, ]
    
    if (nrow(matched_variants) > 0) {
      # Combine information from cleaned_peptides and matched_variants
      for (j in 1:nrow(matched_variants)) {
        result <- cbind(cleaned_peptides[i, ], 
                        matched_variants[j, c("Consequence", "BIOTYPE", "VAF")])
        matched_results[[length(matched_results) + 1]] <- result
      }
    }
  }
}

# Combine all results into a single data frame
if (length(matched_results) > 0) {
  final_result <- do.call(rbind, matched_results)
  
  # Save the final result to a CSV file
  write.csv(final_result, "matched_features_result.csv", row.names = FALSE)
} else {
  cat("No matches found\n")
}



getwd()
setwd ("C:/Users/agsko/dev/pcm/predictions/")

# Load your CSV files
merged_genes <- read.csv("C:/Users/agsko/dev/pcm/predictions/merged_common_genes_upregulated_with_groups.csv")
matched_features <- read.csv("C:/Users/agsko/dev/pcm/predictions/matched_features_result.csv")

# Perform the matching by 'Gene' and keep only matched rows
matched_data <- merge(matched_features, merged_genes[, c("Gene", "DEG_Group")], by = "Gene", all = FALSE)

# Save the matched result to a new CSV file
write.csv(matched_data, "C:/Users/agsko/dev/pcm/predictions/final_matched_genes_with_deg_group_filtered.csv", row.names = FALSE)

# Print the first few rows to check
head(matched_data)




# Load your CSV files
merged_genes <- read.csv("C:/Users/agsko/dev/pcm/predictions/merged_upregulated_genes_unique.csv")
matched_features <- read.csv("C:/Users/agsko/dev/pcm/predictions/matched_features_result.csv")

# Perform the matching by 'Gene' and keep only matched rows
matched_data <- merge(matched_features, merged_genes[, c("Gene", "DEG_Group")], by = "Gene", all = FALSE)

# Save the matched result to a new CSV file
write.csv(matched_data, "C:/Users/agsko/dev/pcm/predictions/final_matched_genes_unique.csv", row.names = FALSE)

# List of files and corresponding identifiers for origin
file_paths <- c(
  "C:/Users/agsko/dev/pcm/EdgeR/Venn and common genes/unique_ure_rt_dmxaa_genes_upregulated_with_names.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/Venn and common genes/unique_ure_rt_genes_upregulated_with_names.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/Venn and common genes/unique_ure_dmxaa_genes_upregulated_with_names.csv",
  "C:/Users/agsko/dev/pcm/EdgeR/Venn and common genes/unique_ure_genes_upregulated_with_names.csv"
)

origin_labels <- c("unique_ure_rt_dmxaa_genes", "unique_ure_rt_genes", "unique_ure_dmxaa_genes", "unique_ure_genes")

# Initialize an empty list to store data frames
combined_data <- list()

# Loop through the files and add the origin column
for (i in 1:length(file_paths)) {
  # Read the CSV file
  df <- read.csv(file_paths[i])
  
  # Add a new column for the origin
  df$Origin <- origin_labels[i]
  
  # Append to the list
  combined_data[[i]] <- df
}

# Combine all data frames into one
final_combined_data <- do.call(rbind, combined_data)

# Save the combined data to a CSV file
write.csv(final_combined_data, "C:/Users/agsko/dev/pcm/predictions/merged_upregulated_genes_unique.csv", row.names = FALSE)

# Print the first few rows to check the output
head(final_combined_data)

#MANUALLY change Origin to DEG_Group


# Define file paths
peptide_file_path <- "C:/Users/agsko/dev/pcm/predictions/combined_class2_predictions_IC500.csv"
combined_data_file_path <- "C:/Users/agsko/dev/pcm/new_lifting/annotated_vcf/parsed_csv_files/combined_data.csv"
output_file_path <- "C:/Users/agsko/dev/pcm/predictions/matched_class2_all_peptides.csv"

# Read the peptides file
peptides_data <- read.csv(peptide_file_path, stringsAsFactors = FALSE)

# Read the combined data with sequences
combined_data <- read.csv(combined_data_file_path, stringsAsFactors = FALSE)

# Initialize an empty data frame to store matched peptides
matched_data <- data.frame()

# Loop through each peptide and search for it as a substring in the combined_data sequences
for (i in 1:nrow(peptides_data)) {
  peptide <- peptides_data$peptide[i]
  
  # Check if the peptide exists as a substring in any sequence in combined_data
  matches <- combined_data[grep(peptide, combined_data$Sequence), ]
  
  # If matches are found, add the relevant information to the matched_data
  if (nrow(matches) > 0) {
    # Add the matched rows along with the corresponding peptide data
    for (j in 1:nrow(matches)) {
      matched_row <- cbind(peptides_data[i, ], matches[j, c("Gene","Transcript", "Type", "Mutation_Type", "Sample")])
      matched_data <- rbind(matched_data, matched_row)
    }
  }
}

# Save the matched data to a new CSV file
write.csv(matched_data, output_file_path, row.names = FALSE)


cat("Matched peptides file saved successfully!\n")


# Define the file path
file_path <- "C:/Users/agsko/dev/pcm/predictions/matched_class2_all_peptides.csv"

# Load the CSV file into a DataFrame
df <- read.csv(file_path)

# Check for duplicates based on all columns
duplicates <- df %>%
  filter(duplicated(.) | duplicated(., fromLast = TRUE))  # Find all duplicate rows

# Display duplicates
if (nrow(duplicates) > 0) {
  print("Found duplicates:")
  print(duplicates)
} else {
  print("No duplicates found.")
}

# Optionally, save the duplicates to a separate file
duplicates_file <- "C:/Users/agsko/dev/pcm/predictions/duplicate_class2_all_peptides.csv"
write.csv(duplicates, file = duplicates_file, row.names = FALSE)
print(paste("Duplicate rows saved to", duplicates_file))

# Define the file path
file_path <- "C:/Users/agsko/dev/pcm/predictions/matched_class2_all_peptides.csv"

# Load the CSV file into a DataFrame
df <- read.csv(file_path)

# Remove the second occurrence of duplicate rows (keep only the first occurrence)
df_clean <- df %>%
  distinct()  # This removes duplicate rows based on all columns

# Save the clean dataset to a new file
clean_file <- "C:/Users/agsko/dev/pcm/predictions/cleaned_class2_all_peptides.csv"
write.csv(df_clean, file = clean_file, row.names = FALSE)

# Load your CSV file
clean_file <- "C:/Users/agsko/dev/pcm/predictions/cleaned_class2_all_peptides.csv"
cleaned_peptides <- read.csv(clean_file)

# Add the Sample_Group column based on the Sample number
cleaned_peptides$Sample_Group <- with(cleaned_peptides, 
                                      ifelse(Sample %in% 1:4, "Ure_mock_IR", 
                                             ifelse(Sample %in% 5:8, "Ure_mock_IR_DMXAA", 
                                                    ifelse(Sample %in% 9:12, "Ure_IR", 
                                                           ifelse(Sample == 13, "PBS_IR_DMXAA", NA)))))

# Save the updated file under the same name
write.csv(cleaned_peptides, clean_file, row.names = FALSE)

# Print the first few rows to verify
head(cleaned_peptides)


# Print the number of rows before and after cleaning
print(paste("Original number of rows:", nrow(df)))
print(paste("Number of rows after removing duplicates:", nrow(df_clean)))

# Print a message
print(paste("Cleaned dataset saved to", clean_file))

# MANUALLY change Transcript to Feature column

library(stringr)


# Load your CSV files
cleaned_peptides <- read.csv("C:/Users/agsko/dev/pcm/predictions/cleaned_class2_all_peptides.csv")
final_variants <- read.csv("C:/Users/agsko/dev/pcm/predictions/final_variants_mm10.csv")

# Function to perform the partial matching
find_partial_matches <- function(peptide_feature, variants_features) {
  # Check if any variant feature contains the peptide feature as a substring
  matches <- grepl(peptide_feature, variants_features)
  return(matches)
}

# Initialize an empty list to store results
matched_results <- list()

# Loop through the cleaned_peptides features
for (i in 1:nrow(cleaned_peptides)) {
  peptide_feature <- cleaned_peptides$Feature[i]
  sample_group <- cleaned_peptides$Sample_Group[i]
  
  # Find matches in final_variants Feature column
  matches <- find_partial_matches(peptide_feature, final_variants$Feature)
  
  # If matches are found
  if (any(matches)) {
    # Subset final_variants based on the matches
    matched_variants <- final_variants[matches & final_variants$Sample_Group == sample_group, ]
    
    if (nrow(matched_variants) > 0) {
      # Combine information from cleaned_peptides and matched_variants
      for (j in 1:nrow(matched_variants)) {
        result <- cbind(cleaned_peptides[i, ], 
                        matched_variants[j, c("Consequence", "BIOTYPE", "VAF")])
        matched_results[[length(matched_results) + 1]] <- result
      }
    }
  }
}

# Combine all results into a single data frame
if (length(matched_results) > 0) {
  final_result <- do.call(rbind, matched_results)
  
  # Save the final result to a CSV file
  write.csv(final_result, "matched_features_class2_result.csv", row.names = FALSE)
} else {
  cat("No matches found\n")
}
getwd()
setwd ("C:/Users/agsko/dev/pcm/predictions/")

# Load your CSV files
merged_genes <- read.csv("C:/Users/agsko/dev/pcm/predictions/merged_common_genes_upregulated_with_groups.csv")
matched_features <- read.csv("C:/Users/agsko/dev/pcm/predictions/matched_features_class2_result.csv")

# Perform the matching by 'Gene' and keep only matched rows
matched_data <- merge(matched_features, merged_genes[, c("Gene", "DEG_Group")], by = "Gene", all = FALSE)

# Save the matched result to a new CSV file
write.csv(matched_data, "C:/Users/agsko/dev/pcm/predictions/final_matched_class2_genes_with_deg_group_filtered.csv", row.names = FALSE)

# Print the first few rows to check
head(matched_data)



#MANUALLY cgane Gene-name to Gene in the merged_uprefulated file
# Load your CSV files
merged_genes <- read.csv("C:/Users/agsko/dev/pcm/predictions/merged_upregulated_genes_unique.csv")
matched_features <- read.csv("C:/Users/agsko/dev/pcm/predictions/matched_features_class2_result.csv")

# Perform the matching by 'Gene' and keep only matched rows
matched_data <- merge(matched_features, merged_genes[, c("Gene", "DEG_Group")], by = "Gene", all = FALSE)

# Save the matched result to a new CSV file
write.csv(matched_data, "C:/Users/agsko/dev/pcm/predictions/final_matched_class2_genes_unique.csv", row.names = FALSE)


########FINAL LISTS

# Define the file paths
file1 <- "C:/Users/agsko/dev/pcm/predictions/final_matched_genes_unique.csv"
file2 <- "C:/Users/agsko/dev/pcm/predictions/final_matched_genes_with_deg_group_filtered.csv"

# Load the CSV files into DataFrames
df1 <- read.csv(file1)
df2 <- read.csv(file2)

# Combine the two DataFrames by rows (assuming they have the same structure)
combined_df <- rbind(df1, df2)

# Save the combined DataFrame to a new CSV file
write.csv(combined_df, "C:/Users/agsko/dev/pcm/predictions/combined_final_matched_genes_class1.csv", row.names = FALSE)

# Print the first few rows to verify
head(combined_df)

# Define the file paths for class 2 genes
file1 <- "C:/Users/agsko/dev/pcm/predictions/final_matched_class2_genes_with_deg_group_filtered.csv"
file2 <- "C:/Users/agsko/dev/pcm/predictions/final_matched_class2_genes_unique.csv"

# Load the CSV files into DataFrames
df1 <- read.csv(file1)
df2 <- read.csv(file2)

# Combine the two DataFrames by rows (assuming they have the same structure)
combined_df <- rbind(df1, df2)

# Save the combined DataFrame to a new CSV file
write.csv(combined_df, "C:/Users/agsko/dev/pcm/predictions/combined_final_matched_class2_genes.csv", row.names = FALSE)

# Print the first few rows to verify
head(combined_df)

# Load the two files
matched_features_file <- "C:/Users/agsko/dev/pcm/predictions/matched_features_result.csv"
top_contributing_genes_file <- "C:/Users/agsko/dev/pcm/predictions/top_contributing_genes_PC1_PC2_80percent_pred.csv"

# Read the CSV files
matched_features_df <- read.csv(matched_features_file)
top_contributing_genes_df <- read.csv(top_contributing_genes_file)

# Filter for matching genes in the 'Gene' column
matched_genes_df <- inner_join(matched_features_df, top_contributing_genes_df, by = "Gene")

# Save the filtered matched genes to a new file
output_file <- "C:/Users/agsko/dev/pcm/predictions/matched_class1_top_contributing_genes.csv"
write.csv(matched_genes_df, output_file, row.names = FALSE)

cat("Filtered matched genes saved to", output_file, "\n")

#Now the same for class 2

# Load the two files
matched_features_file <- "C:/Users/agsko/dev/pcm/predictions/matched_features_class2_result.csv"
top_contributing_genes_file <- "C:/Users/agsko/dev/pcm/predictions/top_contributing_genes_PC1_PC2_80percent_pred.csv"

# Read the CSV files
matched_features_df <- read.csv(matched_features_file)
top_contributing_genes_df <- read.csv(top_contributing_genes_file)

# Filter for matching genes in the 'Gene' column
matched_genes_df <- inner_join(matched_features_df, top_contributing_genes_df, by = "Gene")

# Save the filtered matched genes to a new file
output_file <- "C:/Users/agsko/dev/pcm/predictions/matched_class2_top_contributing_genes.csv"
write.csv(matched_genes_df, output_file, row.names = FALSE)

cat("Filtered matched genes saved to", output_file, "\n")


