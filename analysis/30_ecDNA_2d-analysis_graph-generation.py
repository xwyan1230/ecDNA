import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import shared.dataframe as dat
import numpy as np
import random

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220407_sp8_DMandHSR/"
colors = [(0.8, 0.8, 0.8), (0.85, 0.35, 0.25)]
# colors = [(0.2, 0.2, 0.2), (0.85, 0.35, 0.25)]
rmax = 100
# ['#FFA500', '#40E0D0']

data_DM = pd.read_csv('%sDM.txt' % master_folder, na_values=['.'], sep='\t')
data_HSR = pd.read_csv('%sHSR.txt' % master_folder, na_values=['.'], sep='\t')


# heat map
"""data_DM['sample'] = [1] * len(data_DM)
data_HSR['sample'] = [0] * len(data_HSR)

data = pd.concat([data_DM, data_HSR], axis=0, ignore_index=True)

data_feature = data.copy()
data_feature = data_feature.drop(['FOV', 'nuclear_label', 'nuclear_centroid', 'ecDNA_area', 'ecDNA_mean_int',
                                        'ecDNA_intensity', 'ecDNA_centroid', 'ecDNA_localization_from_centroid',
                                        'ecDNA_distance_from_centroid', 'g', 'dg'], axis=1)

scaler = StandardScaler()
data_feature_scale = pd.DataFrame(scaler.fit_transform(data_feature))
data_feature_scale.columns = data_feature.columns

plt.subplots(figsize=(6, 50))
ax1 = sns.heatmap(data_feature_scale, cbar=0, linewidths=2, vmax=3, vmin=-2, square=True, cmap='viridis')
plt.savefig('%s/heatmap.pdf' % master_folder)
plt.close()"""

# features-number (pyplot probability)
"""features = ['nuclear_area', 'nuclear_major_axis', 'nuclear_minor_axis', 'nuclear_axis_ratio', 'nuclear_circularity',
            'nuclear_eccentricity', 'nuclear_FISH_mean_intensity', 'nuclear_total_intensity', 'ecDNA_number',
            'ecDNA_total_area', 'area_ratio', 'ecDNA_mean_area', 'ecDNA_max_area', 'ecDNA_total_intensity',
            'ecDNA_participating_coefficient', 'MYC_mean_intensity', 'MYC_total_intensity']

for i in features:
    plt.subplots(figsize=(6, 4))

    weights_DM = np.ones_like(data_DM[i]) / len(data_DM)
    weights_HSR = np.ones_like(data_HSR[i]) / len(data_HSR)
    plt.hist([data_HSR[i], data_DM[i]], weights=[weights_HSR, weights_DM], color=colors,
             edgecolor=(0.2, 0.2, 0.2), label=['HSR', 'DM'], bins=20)
    plt.xlabel(i)
    plt.ylabel('Probability')
    plt.legend()
    plt.savefig('%sfeature_%s.pdf' % (master_folder, i))
    plt.close()"""

# features-number (seaborn)
"""data_DM['sample'] = ['DM'] * len(data_DM)
data_HSR['sample'] = ['HSR'] * len(data_HSR)

features = ['nuclear_area', 'nuclear_major_axis', 'nuclear_minor_axis', 'nuclear_axis_ratio', 'nuclear_circularity',
            'nuclear_eccentricity', 'nuclear_FISH_mean_intensity', 'nuclear_total_intensity', 'ecDNA_number',
            'ecDNA_total_area', 'area_ratio', 'ecDNA_mean_area', 'ecDNA_max_area', 'ecDNA_total_intensity',
            'ecDNA_participating_coefficient', 'MYC_mean_intensity', 'MYC_total_intensity']

min_n = min(len(data_DM), len(data_HSR))
lst_DM = random.sample(list(np.arange(len(data_DM))), min_n)
lst_DM.sort()
lst_HSR = random.sample(list(np.arange(len(data_HSR))), min_n)
lst_HSR.sort()
data_DM_select = data_DM.iloc[lst_DM].reset_index()
data_HSR_select = data_HSR.iloc[lst_HSR].reset_index()

data = pd.concat([data_HSR_select, data_DM_select], axis=0, ignore_index=True)

sns.set_palette(sns.color_palette(colors))
for i in features:
    plt.subplots(figsize=(6, 4))
    # ax1 = sns.histplot(data, x='MYC_total_intensity', hue='sample', multiple='dodge', stat='probability', bins=20)
    ax1 = sns.histplot(data, x=i, hue='sample', bins=20, kde=True)
    plt.savefig('%sfeature_seaborn_%s.pdf' % (master_folder, i))
    plt.close()"""

# features-list (pyplot probability)
"""data_DM['sample'] = ['DM'] * len(data_DM)
data_HSR['sample'] = ['HSR'] * len(data_HSR)

features = ['ecDNA_area', 'ecDNA_mean_int', 'ecDNA_intensity', 'ecDNA_distance_from_centroid']

for i in features:
    data_DM[i] = [dat.str_to_float(data_DM[i][j]) for j in range(len(data_DM))]
    data_HSR[i] = [dat.str_to_float(data_HSR[i][j]) for j in range(len(data_HSR))]
    temp = pd.DataFrame()
    feature_DM = []
    feature_HSR = []
    for j in range(len(data_DM)):
        feature_DM = feature_DM + data_DM[i][j]
    for j in range(len(data_HSR)):
        feature_HSR = feature_HSR + data_HSR[i][j]
    temp[i] = feature_DM + feature_HSR
    temp['sample'] = ['DM']*len(feature_DM) + ['HSR']*len(feature_HSR)

    plt.subplots(figsize=(6, 4))
    weights_DM = np.ones_like(feature_DM) / len(feature_DM)
    weights_HSR = np.ones_like(feature_HSR) / len(feature_HSR)
    plt.hist([feature_HSR, feature_DM], weights=[weights_HSR, weights_DM], color=colors,
             edgecolor=(0.2, 0.2, 0.2), label=['HSR', 'DM'], bins=20)
    plt.xlabel(i)
    plt.ylabel('Probability')
    plt.legend()
    plt.savefig('%sfeature_%s.pdf' % (master_folder, i))
    plt.close()"""

# features-list (seaborn)

data_DM['sample'] = ['DM'] * len(data_DM)
data_HSR['sample'] = ['HSR'] * len(data_HSR)

features = ['ecDNA_area', 'ecDNA_mean_int', 'ecDNA_intensity', 'ecDNA_distance_from_centroid']

sns.set_palette(sns.color_palette(colors))
for i in features:
    data_DM[i] = [dat.str_to_float(data_DM[i][j]) for j in range(len(data_DM))]
    data_HSR[i] = [dat.str_to_float(data_HSR[i][j]) for j in range(len(data_HSR))]
    temp = pd.DataFrame()
    feature_DM = []
    feature_HSR = []
    for j in range(len(data_DM)):
        feature_DM = feature_DM + data_DM[i][j]
    for j in range(len(data_HSR)):
        feature_HSR = feature_HSR + data_HSR[i][j]

    min_n = min(len(feature_DM), len(feature_HSR))
    lst_DM = random.sample(feature_DM, min_n)
    lst_HSR = random.sample(feature_HSR, min_n)

    temp[i] = lst_HSR + lst_DM
    temp['sample'] = ['HSR']*len(lst_HSR) + ['DM']*len(lst_DM)

    plt.subplots(figsize=(6, 4))
    # ax1 = sns.histplot(data, x='MYC_total_intensity', hue='sample', multiple='dodge', stat='probability', bins=20)
    ax1 = sns.histplot(temp, x=i, hue='sample', bins=20, kde=True)
    plt.savefig('%sfeature_seaborn_%s.pdf' % (master_folder, i))
    plt.close()

# features comparison
# note, number needs to be similar between sample, otherwise add in random sampling step
"""data_DM['sample'] = ['DM'] * len(data_DM)
data_HSR['sample'] = ['HSR'] * len(data_HSR)

data = pd.concat([data_HSR, data_DM], axis=0, ignore_index=True)

features = ['nuclear_area', 'nuclear_major_axis', 'nuclear_minor_axis', 'nuclear_axis_ratio', 'nuclear_circularity',
            'nuclear_eccentricity', 'nuclear_FISH_mean_intensity', 'nuclear_total_intensity', 'ecDNA_number',
            'ecDNA_total_area', 'area_ratio', 'ecDNA_mean_area', 'ecDNA_max_area', 'ecDNA_total_intensity',
            'ecDNA_participating_coefficient', 'MYC_mean_intensity', 'MYC_total_intensity']
target_feature = 'MYC_mean_intensity'
sns.set_palette(sns.color_palette(colors))
for i in features:
    if i != target_feature:
        # ax1 = sns.histplot(data, x='MYC_total_intensity', y='ecDNA_total_intensity', hue='sample', multiple='dodge',
        # stat='probability', bins=20)
        ax1 = sns.jointplot(data=data, x=target_feature, y=i, hue='sample')
        plt.savefig('%scomparison_of_%s_vs_%s.pdf' % (master_folder, target_feature, i))
        plt.close()"""

# auto-correlation
"""data_DM['g'] = [dat.str_to_float(data_DM['g'][i]) for i in range(len(data_DM))]
data_HSR['g'] = [dat.str_to_float(data_HSR['g'][i]) for i in range(len(data_HSR))]

number_nuclear_DM = len(data_DM)
number_nuclear_HSR = len(data_HSR)

FISH_mean_curve_DM, FISH_ci_lower_DM, FISH_ci_higher_DM = dat.mean_list(data_DM['g'].tolist())
FISH_mean_curve_HSR, FISH_ci_lower_HSR, FISH_ci_higher_HSR = dat.mean_list(data_HSR['g'].tolist())

r = np.arange(0, rmax + 1, 1)

plt.subplots(figsize=(6, 4))
for i in range(len(data_HSR)):
    plt.plot(r, data_HSR['g'][i], alpha=0.05, color=[colors[0][j]+0.1 for j in range(len(colors[0]))])
for i in range(len(data_DM)):
    plt.plot(r, data_DM['g'][i], alpha=0.05, color=[colors[1][j]+0.1 for j in range(len(colors[1]))])
plt.plot(r, FISH_mean_curve_DM, color=colors[1], label='DM, n=%s' % number_nuclear_DM)
plt.plot(r, FISH_ci_lower_DM, color=colors[1], linestyle='--', linewidth=0.5)
plt.plot(r, FISH_ci_higher_DM, color=colors[1], linestyle='--', linewidth=0.5)
plt.plot(r, FISH_mean_curve_HSR, color=colors[0], label='HSR, n=%s' % number_nuclear_HSR)
plt.plot(r, FISH_ci_lower_HSR, color=colors[0], linestyle='--', linewidth=0.5)
plt.plot(r, FISH_ci_higher_HSR, color=colors[0], linestyle='--', linewidth=0.5)
plt.xlabel('r')
plt.ylabel('g')
plt.legend()
plt.ylim([-0.5, 10.5])
plt.savefig('%s/auto_correlation_comparison.pdf' % master_folder)
plt.close()

plt.subplots(figsize=(6, 4))
for i in range(len(data_HSR)):
    plt.plot(r, data_HSR['g'][i], alpha=0.05, color=[colors[0][j]+0.1 for j in range(len(colors[0]))])
plt.plot(r, FISH_mean_curve_HSR, color=colors[0], label='HSR, n=%s' % number_nuclear_HSR)
plt.plot(r, FISH_ci_lower_HSR, color=colors[0], linestyle='--', linewidth=0.5)
plt.plot(r, FISH_ci_higher_HSR, color=colors[0], linestyle='--', linewidth=0.5)
plt.xlabel('r')
plt.ylabel('g')
plt.legend()
plt.ylim([-0.5, 10.5])
plt.savefig('%s/auto_correlation_HSR.pdf' % master_folder)
plt.close()

plt.subplots(figsize=(6, 4))
for i in range(len(data_DM)):
    plt.plot(r, data_DM['g'][i], alpha=0.05, color=[colors[1][j]+0.1 for j in range(len(colors[1]))])
plt.plot(r, FISH_mean_curve_DM, color=colors[1], label='DM, n=%s' % number_nuclear_DM)
plt.plot(r, FISH_ci_lower_DM, color=colors[1], linestyle='--', linewidth=0.5)
plt.plot(r, FISH_ci_higher_DM, color=colors[1], linestyle='--', linewidth=0.5)
plt.xlabel('r')
plt.ylabel('g')
plt.legend()
plt.ylim([-0.5, 10.5])
plt.savefig('%s/auto_correlation_DM.pdf' % master_folder)
plt.close()"""

