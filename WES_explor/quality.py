import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Specify the directory where you want to save the images
output_directory = '/mnt/c/Users/agsko/dev/pcm/WES_explor'

# Ensure the output directory exists
if not os.path.exists(output_directory):
    os.makedirs(output_directory, exist_ok=True)

# Load your CSV file into a pandas DataFrame
df = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/merged_files/final_variants_mm10.csv')

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Specify the directory where you want to save the images
output_directory = '/mnt/c/Users/agsko/dev/pcm/WES_explor'

# Ensure the output directory exists
if not os.path.exists(output_directory):
    os.makedirs(output_directory, exist_ok=True)

# Load your CSV file into a pandas DataFrame
df = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/merged_files/final_variants_mm10.csv')

# Calculate summary statistics for mapping quality columns
mapping_quality_stats = df[['MQ_x', 'MQ_y', 'MQ0']].describe()
print(mapping_quality_stats)

# Filter rows with high MQ0 values (above a certain threshold, e.g., 1000)
high_mq0 = df[df['MQ0'] > 1000]
print(high_mq0)

# Check if the QSS column exists
if 'QSS' in df.columns:
    # Basic statistics for the QSS column
    qss_stats = df['QSS'].describe()

    # Print the results
    print("QSS Statistics:")
    print(qss_stats)
else:
    print("QSS column not found in the dataset.")

# Optional: If you want to save the statistics to a file
qss_stats.to_csv('/mnt/c/Users/agsko/dev/pcm/WES_explor/qss_statistics.csv')

# Calculate the percentage of QSS scores greater than 30
qss_above_30 = df['QSS'][df['QSS'] > 30].count()
total_qss = df['QSS'].count()
percentage_above_30 = (qss_above_30 / total_qss) * 100

# Calculate the median QSS score
median_qss = df['QSS'].median()

# Print the results
print(f"Percentage of QSS scores above 30: {percentage_above_30:.2f}%")
print(f"Median QSS score: {median_qss}")

import pandas as pd
from scipy.stats import chi2_contingency

# Load your CSV file into a pandas DataFrame
df = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/merged_files/final_variants_mm10.csv')

# Set the threshold for "high" MQ0 and QSS values
mq0_threshold = 50  # Adjust this based on your data distribution
qss_threshold = 30  # Common threshold for high-quality variants

# Filter for variants with high MQ0 and high QSS
high_mq0_high_qss_variants = df[(df['MQ0'] > mq0_threshold) & (df['QSS'] > qss_threshold)]

# Output the filtered variants
print("Variants with High MQ0 and High QSS:")
print(high_mq0_high_qss_variants[['CHROM_x', 'POS', 'REF_x', 'ALT_x', 'MQ0', 'QSS', 'IMPACT', 'Consequence']])


# Filter the data for control + PBS and ure-treated samples
control_pbs = df[df['Sample_Group'].isin(['control', 'PBS'])]
ure_treated = df[df['Sample_Group'].str.contains('Ure')]

# Create a contingency table of counts for Variant_Type across the two groups
contingency_table = pd.crosstab(df['Sample_Group'].apply(lambda x: 'Control_PBS' if x in ['control', 'PBS'] else 'Ure_Treated'),
                                df['Variant_Type'])

# Perform the chi-square test
chi2, p, dof, expected = chi2_contingency(contingency_table)

# Output the results
print("Chi-square statistic:", chi2)
print("p-value:", p)
print("Degrees of freedom:", dof)
print("Expected frequencies:", expected)

# Interpret the result
if p < 0.05:
    print("There is a statistically significant difference in the counts of SNVs and indels between control + PBS and ure-treated samples.")
else:
    print("There is no statistically significant difference in the counts of SNVs and indels between control + PBS and ure-treated samples.")
