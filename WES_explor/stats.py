import pandas as pd
import scipy.stats as stats  # Ensure that scipy.stats is imported
import matplotlib.pyplot as plt
import seaborn as sns

# Load the frameshift summary data
frameshift_summary = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/WES_explor/frameshift_summary_by_group.csv')

# Perform Kruskal-Wallis Test
kruskal_stat, p_value = stats.kruskal(*[frameshift_summary[frameshift_summary['Sample_Group'] == group]['Frameshift_Indels_Count'] for group in frameshift_summary['Sample_Group'].unique()])
print(f"Kruskal-Wallis Test: H-statistic={kruskal_stat}, p-value={p_value}")

if p_value < 0.05:
    print("There is a statistically significant difference in frameshift indels counts among the groups.")
else:
    print("No statistically significant difference in frameshift indels counts among the groups.")


import pandas as pd
import scipy.stats as stats

# Load the summary data
frameshift_summary = pd.read_csv('/mnt/c/Users/agsko/dev/pcm/WES_explor/frameshift_summary_by_group.csv')

# Combine the Urethane-treated groups into a single group
ure_groups = frameshift_summary[frameshift_summary['Sample_Group'].isin(['Ure_IR', 'Ure_mock_IR', 'Ure_mock_IR_DMXAA'])]['Frameshift_Indels_Count']
pbs_group = frameshift_summary[frameshift_summary['Sample_Group'] == 'PBS_IR_DMXAA']['Frameshift_Indels_Count']

# Perform the Mann-Whitney U Test
mannwhitney_stat, p_value = stats.mannwhitneyu(ure_groups, pbs_group)

# Output the results
print(f"Mann-Whitney U Test: U-statistic={mannwhitney_stat}, p-value={p_value}")
if p_value < 0.05:
    print("Statistically significant difference in frameshift indels counts between Urethane-treated groups and PBS group.")
else:
    print("No statistically significant difference in frameshift indels counts between Urethane-treated groups and PBS group.")
