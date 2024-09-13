# Load necessary library
library(glue)

# Function to parse mutant peptide files correctly
parse_mutant_peptides_v3 <- function(file_path, output_file) {
  # Read the mutant peptides file
  mutant_peptides <- readLines(file_path)
  
  # Initialize empty lists to store parsed information
  gene_names <- c()
  transcript_ids <- c()
  peptides <- c()
  wt_mt <- c()
  
  # Temporary variables to store current peptide information
  current_peptide <- ""
  current_gene_name <- ""
  current_transcript_id <- ""
  current_wt_mt <- ""
  
  for (line in mutant_peptides) {
    # Check if the line is a header (starts with '>')
    if (grepl("^>", line)) {
      # If we already have a peptide from previous lines, save it
      if (current_peptide != "") {
        # Store the previous peptide information
        gene_names <- c(gene_names, current_gene_name)
        transcript_ids <- c(transcript_ids, current_transcript_id)
        peptides <- c(peptides, current_peptide)
        wt_mt <- c(wt_mt, current_wt_mt)
        
        # Reset current peptide
        current_peptide <- ""
      }
      
      # Extract gene name, transcript ID, and WT/MT status from the header line
      current_gene_name <- gsub(".*\\.([^.]+)\\.ENSMUST[0-9]+.*", "\\1", line)
      current_transcript_id <- gsub(".*\\.(ENSMUST[0-9]+)\\..*", "\\1", line)
      current_wt_mt <- ifelse(grepl("^>WT", line), "WT", "MT")
      
    } else {
      # Concatenate the sequence lines to form the full peptide sequence
      current_peptide <- paste0(current_peptide, line)
    }
  }
  
  # Handle the last peptide in the file (if exists)
  if (current_peptide != "") {
    gene_names <- c(gene_names, current_gene_name)
    transcript_ids <- c(transcript_ids, current_transcript_id)
    peptides <- c(peptides, current_peptide)
    wt_mt <- c(wt_mt, current_wt_mt)
  }
  
  # Create a data frame with the extracted information
  mutant_peptides_df <- data.frame(
    Gene.name = gene_names,
    Transcript.ID = transcript_ids,
    Peptide = peptides,
    WT_MT = wt_mt,
    stringsAsFactors = FALSE
  )
  
  # Save the parsed data to a CSV file for review
  write.csv(mutant_peptides_df, output_file, row.names = FALSE)
  
  # Print a message indicating that the file has been created
  cat(glue("Parsed mutant peptides saved as '{output_file}'.\n"))
}

# Example usage for S4 file
file_path_s4 <- "C:/Users/agsko/dev/pcm/samples_peptides/Sample_2_S4/Sample_2_S4_indels_output_class2_mutant_peptides.txt"
output_file_s4 <- "parsed_mutant_peptides_S4_refined_v3.csv"
parse_mutant_peptides_v3(file_path_s4, output_file_s4)

# Define file paths for all samples
file_paths <- c(
  "C:/Users/agsko/dev/pcm/samples_peptides/Sample_1_S3/Sample_1_S3_indels_output_class2_mutant_peptides.txt",
  "C:/Users/agsko/dev/pcm/samples_peptides/Sample_2_S4/Sample_2_S4_indels_output_class2_mutant_peptides.txt",
  "C:/Users/agsko/dev/pcm/samples_peptides/Sample_3_S5/Sample_3_S5_indels_output_class2_mutant_peptides.txt",
  "C:/Users/agsko/dev/pcm/samples_peptides/Sample_4_S6/Sample_4_S6_indels_output_class2_mutant_peptides.txt",
  "C:/Users/agsko/dev/pcm/samples_peptides/Sample_5_S7/Sample_5_S7_indels_output_class2_mutant_peptides.txt",
  "C:/Users/agsko/dev/pcm/samples_peptides/Sample_6_S8/Sample_6_S8_indels_output_class2_mutant_peptides.txt",
  "C:/Users/agsko/dev/pcm/samples_peptides/Sample_7_S9/Sample_7_S9_indels_output_class2_mutant_peptides.txt",
  "C:/Users/agsko/dev/pcm/samples_peptides/Sample_8_S10/Sample_8_S10_indels_output_class2_mutant_peptides.txt",
  "C:/Users/agsko/dev/pcm/samples_peptides/Sample_9_S11/Sample_9_S11_indels_output_class2_mutant_peptides.txt",
  "C:/Users/agsko/dev/pcm/samples_peptides/Sample_10_S12/Sample_10_S12_indels_output_class2_mutant_peptides.txt",
  "C:/Users/agsko/dev/pcm/samples_peptides/Sample_11_S13/Sample_11_S13_indels_output_class2_mutant_peptides.txt",
  "C:/Users/agsko/dev/pcm/samples_peptides/Sample_12_S14/Sample_12_S14_indels_output_class2_mutant_peptides.txt",
  "C:/Users/agsko/dev/pcm/samples_peptides/Sample_13_S15/Sample_13_S15_indels_output_class2_mutant_peptides.txt"
)

# Define output file names for the parsed files
output_files <- c(
  "parsed_mutant_peptides_S3_refined_v3.csv",
  "parsed_mutant_peptides_S4_refined_v3.csv",
  "parsed_mutant_peptides_S5_refined_v3.csv",
  "parsed_mutant_peptides_S6_refined_v3.csv",
  "parsed_mutant_peptides_S7_refined_v3.csv",
  "parsed_mutant_peptides_S8_refined_v3.csv",
  "parsed_mutant_peptides_S9_refined_v3.csv",
  "parsed_mutant_peptides_S10_refined_v3.csv",
  "parsed_mutant_peptides_S11_refined_v3.csv",
  "parsed_mutant_peptides_S12_refined_v3.csv",
  "parsed_mutant_peptides_S13_refined_v3.csv",
  "parsed_mutant_peptides_S14_refined_v3.csv",
  "parsed_mutant_peptides_S15_refined_v3.csv"
)

# Apply the parsing function to each file and save the results
for (i in seq_along(file_paths)) {
  parse_mutant_peptides_v3(file_paths[i], output_files[i])
}

# Load RNA-seq TPM data
tpm_data <- read.csv("C:/Users/agsko/dev/pcm/TPM_values.csv")

# Filter TPM data for Sample X88 (as an example)
tpm_x88 <- tpm_data %>%
  dplyr::select(Gene.name, X88) %>%
  dplyr::filter(X88 > 0)  # Keep only expressed genes (TPM > 0)


# Define the file path to your mutant peptides file for Sample X88
file_path_x88 <- "C:/Users/agsko/dev/pcm/parsed_mutant_peptides_S4_refined_v3.csv"

# Read the mutant peptides CSV file
mutant_peptides_df <- read.csv(file_path_x88, stringsAsFactors = FALSE)

# Filter the mutant peptides based on the gene expression in X88
# Keep only mutant peptides (WT_MT == "MT") for expressed genes
expressed_peptides <- mutant_peptides_df %>%
  filter(Gene.name %in% tpm_x88$Gene.name & WT_MT == "MT")

# Preview the filtered peptides
head(expressed_peptides)

# Save the filtered peptides to a CSV file if needed
write.csv(expressed_peptides, "filtered_peptides_X88_class2.csv", row.names = FALSE)

setwd("C:/Users/agsko/dev/pcm")

# Load RNA-seq TPM data
tpm_data <- read.csv("C:/Users/agsko/dev/pcm/TPM_values.csv")

# Function to filter mutant peptides based on expressed genes
filter_expressed_peptides <- function(sample_id, mutant_file_path, output_file_name) {
  # Filter TPM data for the specific sample
  tpm_sample <- tpm_data %>%
    dplyr::select(Gene.name, !!sym(sample_id)) %>%
    dplyr::filter(!!sym(sample_id) > 0)  # Keep only expressed genes (TPM > 0)
  
  # Read the mutant peptides CSV file
  mutant_peptides_df <- read.csv(mutant_file_path, stringsAsFactors = FALSE)
  
  # Filter the mutant peptides based on the gene expression in the sample
  expressed_peptides <- mutant_peptides_df %>%
    filter(Gene.name %in% tpm_sample$Gene.name & WT_MT == "MT")
  
  # Save the filtered peptides to a CSV file
  write.csv(expressed_peptides, output_file_name, row.names = FALSE)
}

# Apply the function to each sample
filter_expressed_peptides("X89", "C:/Users/agsko/dev/pcm/parsed_mutant_peptides_S5_refined_v3.csv", "filtered_peptides_X89_class2.csv")
filter_expressed_peptides("X90", "C:/Users/agsko/dev/pcm/parsed_mutant_peptides_S6_refined_v3.csv", "filtered_peptides_X90_class2.csv")
filter_expressed_peptides("X92", "C:/Users/agsko/dev/pcm/parsed_mutant_peptides_S7_refined_v3.csv", "filtered_peptides_X92_class2.csv")
filter_expressed_peptides("X93", "C:/Users/agsko/dev/pcm/parsed_mutant_peptides_S8_refined_v3.csv", "filtered_peptides_X93_class2.csv")
filter_expressed_peptides("X94", "C:/Users/agsko/dev/pcm/parsed_mutant_peptides_S9_refined_v3.csv", "filtered_peptides_X94_class2.csv")
filter_expressed_peptides("X97", "C:/Users/agsko/dev/pcm/parsed_mutant_peptides_S11_refined_v3.csv", "filtered_peptides_X97_class2.csv")
filter_expressed_peptides("X99", "C:/Users/agsko/dev/pcm/parsed_mutant_peptides_S13_refined_v3.csv", "filtered_peptides_X99_class2.csv")
filter_expressed_peptides("X101", "C:/Users/agsko/dev/pcm/parsed_mutant_peptides_S14_refined_v3.csv", "filtered_peptides_X101_class2.csv")


