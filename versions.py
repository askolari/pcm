import pandas as pd

# Load your CSV file into a pandas DataFrame
df = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/combined_variants_with_sample_info.csv')  # Adjust path if needed

# Summarize the count of each variant type within each sample group
variant_summary = df.groupby(['Sample_Group', 'Variant_Type']).size().unstack(fill_value=0)

print(variant_summary)

import scipy.stats as stats

# Use the observed counts to calculate expected frequencies
chi2, p, dof, expected = stats.chi2_contingency(variant_summary)

# Display the expected frequencies
print("Expected Frequencies:")
print(expected)

# Perform the Chi-Square test
chi2_stat, p_val, dof, expected = stats.chi2_contingency(variant_summary)

# Print the results
print(f"Chi-Square Statistic: {chi2_stat}")
print(f"P-value: {p_val}")
print(f"Degrees of Freedom: {dof}")

import pandas as pd
from scipy.stats import chi2_contingency

# Now, focus only on urethane-treated groups
urethane_groups_data = {
    'Variant_Type': ['indel', 'snv'],
    'Ure_IR': [variant_summary.loc['Ure_IR']['indel'], variant_summary.loc['Ure_IR']['snv']],
    'Ure_mock_IR': [variant_summary.loc['Ure_mock_IR']['indel'], variant_summary.loc['Ure_mock_IR']['snv']],
    'Ure_mock_IR_DMXAA': [variant_summary.loc['Ure_mock_IR_DMXAA']['indel'], variant_summary.loc['Ure_mock_IR_DMXAA']['snv']]
}

df_urethane = pd.DataFrame(urethane_groups_data)
df_urethane.set_index('Variant_Type', inplace=True)

# Transpose for chi-square input
urethane_groups = df_urethane.T

# Run the Chi-Square test on the subsetted urethane-treated groups
chi2_stat_urethane, p_val_urethane, dof_urethane, expected_urethane = chi2_contingency(urethane_groups)

# Print the results for urethane-treated groups
print(f"Chi-Square Statistic (Urethane Groups): {chi2_stat_urethane}")
print(f"P-value (Urethane Groups): {p_val_urethane}")
print(f"Degrees of Freedom (Urethane Groups): {dof_urethane}")

# Optionally, print expected frequencies for urethane groups
print("Expected Frequencies (Urethane Groups):")
print(expected_urethane)

import pandas as pd
from scipy.stats import chi2_contingency

# Data for each pairwise comparison
pairwise_data = {
    'Ure_IR_vs_Ure_mock_IR': {
        'data': {'Variant_Type': ['indel', 'snv'], 'Ure_IR': [957, 3821], 'Ure_mock_IR': [953, 4314]},
    },
    'Ure_IR_vs_Ure_mock_IR_DMXAA': {
        'data': {'Variant_Type': ['indel', 'snv'], 'Ure_IR': [957, 3821], 'Ure_mock_IR_DMXAA': [816, 4581]},
    },
    'Ure_mock_IR_vs_Ure_mock_IR_DMXAA': {
        'data': {'Variant_Type': ['indel', 'snv'], 'Ure_mock_IR': [953, 4314], 'Ure_mock_IR_DMXAA': [816, 4581]},
    }
}

# Perform the pairwise Chi-Square tests
pairwise_results = {}
for comparison, details in pairwise_data.items():
    df = pd.DataFrame(details['data'])
    df.set_index('Variant_Type', inplace=True)
    chi2_stat, p_val, dof, expected = chi2_contingency(df.T)
    pairwise_results[comparison] = {'Chi2_Stat': chi2_stat, 'P-Value': p_val, 'DOF': dof, 'Expected Frequencies': expected}

# Print the results
for comparison, results in pairwise_results.items():
    print(f"{comparison}:")
    print(f"  Chi-Square Statistic: {results['Chi2_Stat']}")
    print(f"  P-value: {results['P-Value']}")
    print(f"  Degrees of Freedom: {results['DOF']}")
    print(f"  Expected Frequencies:\n{results['Expected Frequencies']}")
    print()

import matplotlib.pyplot as plt
import seaborn as sns

# Create the bar plots for each pairwise comparison
plt.figure(figsize=(14, 10))

# Plot for Ure_IR vs. Ure_mock_IR
plt.subplot(3, 1, 1)
df_1 = pd.DataFrame(pairwise_data['Ure_IR_vs_Ure_mock_IR']['data'])
df_1.set_index('Variant_Type', inplace=True)
df_1.plot(kind='bar', ax=plt.gca())
plt.title('Comparison of Variant Types: Ure_IR vs. Ure_mock_IR')
plt.ylabel('Count')

# Plot for Ure_IR vs. Ure_mock_IR_DMXAA
plt.subplot(3, 1, 2)
df_2 = pd.DataFrame(pairwise_data['Ure_IR_vs_Ure_mock_IR_DMXAA']['data'])
df_2.set_index('Variant_Type', inplace=True)
df_2.plot(kind='bar', ax=plt.gca())
plt.title('Comparison of Variant Types: Ure_IR vs. Ure_mock_IR_DMXAA')
plt.ylabel('Count')

# Plot for Ure_mock_IR vs. Ure_mock_IR_DMXAA
plt.subplot(3, 1, 3)
df_3 = pd.DataFrame(pairwise_data['Ure_mock_IR_vs_Ure_mock_IR_DMXAA']['data'])
df_3.set_index('Variant_Type', inplace=True)
df_3.plot(kind='bar', ax=plt.gca())
plt.title('Comparison of Variant Types: Ure_mock_IR vs. Ure_mock_IR_DMXAA')
plt.ylabel('Count')

plt.tight_layout()
plt.show()
