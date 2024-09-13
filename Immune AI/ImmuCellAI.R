
# Load necessary libraries
library(readr)

# Step 1: Load the TPM_values_EPIC.txt file (assuming it's tab-delimited)
tpm_values <- read.table("C:/Users/agsko/dev/pcm/TPM_values_EPIC.txt", header = TRUE, sep = "\t")

# Load necessary libraries
library(readr)

# Step 1: Load the TPM_values_EPIC.txt file (assuming it's tab-delimited)
tpm_values <- read.table("TPM_values_EPIC.txt", header = TRUE, sep = "\t")

# Step 2: Exclude the "Gene_Symbol" column as it is not a sample
tpm_values_data <- tpm_values[, -1]  # Remove the 'Gene_Symbol' column

# Step 3: Add the Group information
# Define the group assignment
Group <- c("Ure", "Ure", "Ure", 
           "Ure_DMXAA", "Ure_DMXAA", "Ure_DMXAA", 
           "Ure_RT", "Ure_RT", "Ure_RT", 
           "Ure_RT_DMXAA", "Ure_RT_DMXAA", "Ure_RT_DMXAA", 
           "PBS_RT_DMXAA", "PBS_RT_DMXAA", "PBS_RT_DMXAA")

# Step 4: Create a new data frame for the group row
group_row <- data.frame(t(Group))

# Set the correct column names for the group row based on the sample columns
colnames(group_row) <- colnames(tpm_values_data)

# Step 5: Insert the group row as the first row in the data frame
tpm_values_final <- rbind(group_row, tpm_values_data)

# Step 6: Add the placeholder for the "Gene_Symbol" column for the group row
gene_symbols <- c("Group", tpm_values$Gene_Symbol)

# Step 7: Add back the 'Gene_Symbol' column to the final dataset
tpm_values_final <- cbind(Gene_Symbol = gene_symbols, tpm_values_final)

# Step 8: Save the updated dataset as a new tab-delimited file
write.table(tpm_values_final, "TPM_values_EPIC_with_groups.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Step 9: Check the final dataset
head(tpm_values_final)
