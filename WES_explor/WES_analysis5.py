import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import os

# Define the directory containing the merged files and where to save the output
merged_files_dir = '/mnt/c/Users/agsko/dev/pcm/merged_files/'
output_directory = '/mnt/c/Users/agsko/dev/pcm/WES_explor/common_gene_analysis'

# Ensure the output directory exists
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Load the final dataset into a pandas DataFrame
df = pd.read_csv(os.path.join(merged_files_dir, 'final_variants_mm10.csv'))

# Ensure the SYMBOL column is treated as a string and handle NaN values
df['SYMBOL'] = df['SYMBOL'].astype(str).fillna('')

# Extract unique genes (SYMBOL) for each Ure condition
ure_ir_genes = set(df[df['Sample_Group'] == 'Ure_IR']['SYMBOL'].unique())
ure_mock_ir_genes = set(df[df['Sample_Group'] == 'Ure_mock_IR']['SYMBOL'].unique())
ure_mock_ir_dmxaa_genes = set(df[df['Sample_Group'] == 'Ure_mock_IR_DMXAA']['SYMBOL'].unique())

# Remove any empty strings if present
ure_ir_genes = {gene for gene in ure_ir_genes if gene.strip()}  # Using .strip() to remove any whitespace
ure_mock_ir_genes = {gene for gene in ure_mock_ir_genes if gene.strip()}
ure_mock_ir_dmxaa_genes = {gene for gene in ure_mock_ir_dmxaa_genes if gene.strip()}

# 1. Common genes between Ure_IR and Ure_mock_IR
common_ure_ir_ure_mock_ir = ure_ir_genes.intersection(ure_mock_ir_genes)
common_ure_ir_ure_mock_ir_file = os.path.join(output_directory, 'common_ure_ir_ure_mock_ir_genes.csv')
pd.DataFrame({'Common_Genes': list(common_ure_ir_ure_mock_ir)}).to_csv(common_ure_ir_ure_mock_ir_file, index=False)
print(f"Common genes between Ure_IR and Ure_mock_IR saved to {common_ure_ir_ure_mock_ir_file}")

# 2. Common genes between Ure_IR and Ure_mock_IR_DMXAA
common_ure_ir_ure_mock_ir_dmxaa = ure_ir_genes.intersection(ure_mock_ir_dmxaa_genes)
common_ure_ir_ure_mock_ir_dmxaa_file = os.path.join(output_directory, 'common_ure_ir_ure_mock_ir_dmxaa_genes.csv')
pd.DataFrame({'Common_Genes': list(common_ure_ir_ure_mock_ir_dmxaa)}).to_csv(common_ure_ir_ure_mock_ir_dmxaa_file, index=False)
print(f"Common genes between Ure_IR and Ure_mock_IR_DMXAA saved to {common_ure_ir_ure_mock_ir_dmxaa_file}")

# 3. Common genes between Ure_mock_IR and Ure_mock_IR_DMXAA
common_ure_mock_ir_ure_mock_ir_dmxaa = ure_mock_ir_genes.intersection(ure_mock_ir_dmxaa_genes)
common_ure_mock_ir_ure_mock_ir_dmxaa_file = os.path.join(output_directory, 'common_ure_mock_ir_ure_mock_ir_dmxaa_genes.csv')
pd.DataFrame({'Common_Genes': list(common_ure_mock_ir_ure_mock_ir_dmxaa)}).to_csv(common_ure_mock_ir_ure_mock_ir_dmxaa_file, index=False)
print(f"Common genes between Ure_mock_IR and Ure_mock_IR_DMXAA saved to {common_ure_mock_ir_ure_mock_ir_dmxaa_file}")

# 4. Common genes among all three Ure conditions
common_all_three_ure_conditions = ure_ir_genes.intersection(ure_mock_ir_genes).intersection(ure_mock_ir_dmxaa_genes)
common_all_three_ure_conditions_file = os.path.join(output_directory, 'common_all_three_ure_conditions_genes.csv')
pd.DataFrame({'Common_Genes': list(common_all_three_ure_conditions)}).to_csv(common_all_three_ure_conditions_file, index=False)
print(f"Common genes among all three Ure conditions saved to {common_all_three_ure_conditions_file}")

# Create a Venn diagram comparing the three Ure conditions
plt.figure(figsize=(10, 10))
venn3([ure_ir_genes, ure_mock_ir_genes, ure_mock_ir_dmxaa_genes],
      set_labels=('Ure_IR', 'Ure_mock_IR', 'Ure_mock_IR_DMXAA'),
      set_colors=('gold', 'lightcoral', 'lightblue'))
plt.title("Venn Diagram: Comparing Ure_IR, Ure_mock_IR, and Ure_mock_IR_DMXAA")
plt.savefig(os.path.join(output_directory, 'venn_ure_conditions_comparison.png'))
plt.show()
