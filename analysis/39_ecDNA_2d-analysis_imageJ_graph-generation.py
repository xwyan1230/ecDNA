import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import shared.dataframe as dat
import numpy as np
import random
import shared.math as mat

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220301_ecDNA_ctrlAndJQ1_NatashaFile/"
colors = [(0.2, 0.2, 0.2), (0.85, 0.35, 0.25)]
rmax = 150
radial_interval = 1
radial_max = 120
relative_radial_interval = 0.01

data_sample = pd.read_csv('%sJQ13hr.txt' % master_folder, na_values=['.'], sep='\t')
data_ctrl = pd.read_csv('%sDMSO.txt' % master_folder, na_values=['.'], sep='\t')

# auto-correlation
feature = 'g_correct'
sample = 'JQ13hr'
ctrl = 'DMSO'

data_sample[feature] = [dat.str_to_float(data_sample[feature][i]) for i in range(len(data_sample))]
data_ctrl[feature] = [dat.str_to_float(data_ctrl[feature][i]) for i in range(len(data_ctrl))]

number_nuclear_sample = len(data_sample)
number_nuclear_ctrl = len(data_ctrl)

mean_curve_sample, ci_lower_sample, ci_higher_sample = dat.mean_list(data_sample[feature].tolist())
mean_curve_ctrl, ci_lower_ctrl, ci_higher_ctrl = dat.mean_list(data_ctrl[feature].tolist())

r = np.arange(0, rmax + 1, 1)

plt.subplots(figsize=(6, 4))
for i in range(len(data_ctrl)):
    plt.plot(r, data_ctrl[feature][i], alpha=0.05, color=[colors[0][j]+0.1 for j in range(len(colors[0]))])
for i in range(len(data_sample)):
    plt.plot(r, data_sample[feature][i], alpha=0.05, color=[colors[1][j]+0.1 for j in range(len(colors[1]))])
plt.plot(r, mean_curve_sample, color=colors[1], label='%s, n=%s' % (sample, number_nuclear_sample))
plt.plot(r, ci_lower_sample, color=colors[1], linestyle='--', linewidth=0.5)
plt.plot(r, ci_higher_sample, color=colors[1], linestyle='--', linewidth=0.5)
plt.plot(r, mean_curve_ctrl, color=colors[0], label='%s, n=%s' % (ctrl, number_nuclear_ctrl))
plt.plot(r, ci_lower_ctrl, color=colors[0], linestyle='--', linewidth=0.5)
plt.plot(r, ci_higher_ctrl, color=colors[0], linestyle='--', linewidth=0.5)
plt.xlabel('r')
plt.ylabel(feature)
plt.legend()
plt.ylim([-0.5, 35.5])
plt.savefig('%s/%s_comparison.pdf' % (master_folder, feature))
plt.close()

# vs MYC expression
"""feature = 'g_value'
sample = 'DM'
ctrl = 'HSR'

data_sample['sample'] = [sample] * len(data_sample)
data_ctrl['sample'] = [ctrl] * len(data_ctrl)

data = pd.concat([data_ctrl, data_sample], axis=0, ignore_index=True)
data['MYC_expression_probability_per_copy'] = \
    data['total_intensity_MYC_IF']/data['total_intensity_MYC_DNAFISH_in_nucleus']

sns.set_palette(sns.color_palette(colors))

ax1 = sns.jointplot(data=data, x=feature, y='MYC_expression_probability_per_copy', hue='sample', alpha=0.7, s=10)
plt.savefig('%scomparison_of_%s_vs_MYC_expression_probability_per_copy.pdf' % (master_folder, feature))
plt.close()
ax1 = sns.jointplot(data=data, x=feature, y='total_intensity_MYC_IF', hue='sample', alpha=0.7, s=10)
plt.savefig('%scomparison_of_%s_vs_MYC_expression.pdf' % (master_folder, feature))
plt.close()
ax1 = sns.jointplot(data=data, x=feature, y='total_intensity_MYC_DNAFISH_in_nucleus', hue='sample', alpha=0.7, s=10)
plt.savefig('%scomparison_of_%s_vs_total_intensity_MYC_DNAFISH_in_nucleus.pdf' % (master_folder, feature))
plt.close()"""

