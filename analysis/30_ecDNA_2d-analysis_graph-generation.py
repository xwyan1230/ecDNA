import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import shared.dataframe as dat
import numpy as np
import random
import shared.math as mat

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220407_sp8_DMandHSR/"
# colors = [(0.8, 0.8, 0.8), (0.85, 0.35, 0.25)]
colors = [(0.2, 0.2, 0.2), (0.85, 0.35, 0.25)]
rmax = 100
radial_interval = 1
radial_max = 120
relative_radial_interval = 0.01
# ['#FFA500', '#40E0D0']

data_DM = pd.read_csv('%sDM.txt' % master_folder, na_values=['.'], sep='\t')
data_HSR = pd.read_csv('%sHSR.txt' % master_folder, na_values=['.'], sep='\t')
data_auto = pd.read_csv('%sDM_auto_correlation_fov0_i2.txt' % master_folder, na_values=['.'], sep='\t')

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

"""data_DM['sample'] = ['DM'] * len(data_DM)
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
    plt.close()"""

# features comparison
# note, number needs to be similar between sample, otherwise add in random sampling step
"""data_DM['sample'] = ['DM'] * len(data_DM)
data_HSR['sample'] = ['HSR'] * len(data_HSR)

data = pd.concat([data_HSR, data_DM], axis=0, ignore_index=True)

features = ['nuclear_area', 'nuclear_major_axis', 'nuclear_minor_axis', 'nuclear_axis_ratio', 'nuclear_circularity',
            'nuclear_eccentricity', 'nuclear_FISH_mean_intensity', 'nuclear_total_intensity', 'ecDNA_number',
            'ecDNA_total_area', 'area_ratio', 'ecDNA_mean_area', 'ecDNA_max_area', 'ecDNA_total_intensity',
            'ecDNA_participating_coefficient', 'MYC_mean_intensity', 'MYC_total_intensity']
target_feature = 'MYC_total_intensity'
sns.set_palette(sns.color_palette(colors))
for i in features:
    if i != target_feature:
        # ax1 = sns.histplot(data, x='MYC_total_intensity', y='ecDNA_total_intensity', hue='sample', multiple='dodge',
        # stat='probability', bins=20)
        ax1 = sns.jointplot(data=data, x=target_feature, y=i, hue='sample')
        plt.savefig('%scomparison_of_%s_vs_%s.pdf' % (master_folder, target_feature, i))
        plt.close()"""

# MYC_total_intensity vs FISH_total_intensity_nuclear
"""data_DM['sample'] = ['DM'] * len(data_DM)
data_HSR['sample'] = ['HSR'] * len(data_HSR)

data = pd.concat([data_HSR, data_DM], axis=0, ignore_index=True)
sns.set_palette(sns.color_palette(colors))

ax1 = sns.jointplot(data=data, x='MYC_total_intensity', y='nuclear_total_intensity', hue='sample')

_, int_fit_r2_DM, int_fit_a_DM = \
    mat.fitting_linear_b0(data[data['sample'] == 'DM']['MYC_total_intensity'].tolist(),
                          data[data['sample'] == 'DM']['nuclear_total_intensity'].tolist())
_, int_fit_r2_HSR, int_fit_a_HSR = \
    mat.fitting_linear_b0(data[data['sample'] == 'HSR']['MYC_total_intensity'].tolist(),
                          data[data['sample'] == 'HSR']['nuclear_total_intensity'].tolist())

x = np.arange(0, 1.5 * pow(10, 9), pow(10, 7))
y_DM = int_fit_a_DM * x
y_HSR = int_fit_a_HSR * x
ax1.ax_joint.plot(x, y_HSR, linewidth=2, color=colors[0], linestyle='--')
ax1.ax_joint.plot(x, y_DM, linewidth=2, color=colors[1], linestyle='--')

plt.savefig('%scomparison_of_MYC_total_intensity_vs_FISH_total_intensity_nuclear.pdf' % master_folder)
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

# g_value vs intensity
"""data_DM['sample'] = ['DM'] * len(data_DM)
data_HSR['sample'] = ['HSR'] * len(data_HSR)

data_DM['g'] = [dat.str_to_float(data_DM['g'][i]) for i in range(len(data_DM))]
data_HSR['g'] = [dat.str_to_float(data_HSR['g'][i]) for i in range(len(data_HSR))]

data = pd.concat([data_HSR, data_DM], axis=0, ignore_index=True)
data['g_value'] = [(data['g'][i][1] + data['g'][i][2] + data['g'][i][3] + data['g'][i][4] + data['g'][i][5])*0.2
                   for i in range(len(data))]
data['MYC_intensity/FISH_intensity'] = data['MYC_total_intensity']/data['nuclear_total_intensity']
# data.to_csv('%sprocessed.txt' % master_folder, index=False, sep='\t')

sns.set_palette(sns.color_palette(colors))

ax1 = sns.jointplot(data=data, x='g_value', y='MYC_intensity/FISH_intensity', hue='sample')

_, g_fit_r2_DM, g_fit_a_DM = \
    mat.fitting_linear_b0(data[data['sample'] == 'DM']['g_value'].tolist(),
                          data[data['sample'] == 'DM']['MYC_intensity/FISH_intensity'].tolist())
_, g_fit_r2_HSR, g_fit_a_HSR = \
    mat.fitting_linear_b0(data[data['sample'] == 'HSR']['g_value'].tolist(),
                          data[data['sample'] == 'HSR']['MYC_intensity/FISH_intensity'].tolist())

x = np.arange(0, 20, 0.5)
y_DM = g_fit_a_DM * x
y_HSR = g_fit_a_HSR * x
ax1.ax_joint.plot(x, y_HSR, linewidth=2, color=colors[0], linestyle='--')
ax1.ax_joint.plot(x, y_DM, linewidth=2, color=colors[1], linestyle='--')

plt.savefig('%scomparison_of_g_value_vs_intensity_residue.pdf' % master_folder)
plt.close()"""


"""features = ['nuclear_area', 'nuclear_major_axis', 'nuclear_minor_axis', 'nuclear_axis_ratio', 'nuclear_circularity',
            'nuclear_eccentricity', 'nuclear_FISH_mean_intensity', 'nuclear_total_intensity', 'ecDNA_number',
            'ecDNA_total_area', 'area_ratio', 'ecDNA_mean_area', 'ecDNA_max_area', 'ecDNA_total_intensity',
            'ecDNA_participating_coefficient', 'MYC_mean_intensity', 'MYC_total_intensity']
target_feature = 'g_value'
sns.set_palette(sns.color_palette(colors))
for i in features:
    if i != target_feature:
        # ax1 = sns.histplot(data, x='MYC_total_intensity', y='ecDNA_total_intensity', hue='sample', multiple='dodge',
        # stat='probability', bins=20)
        ax1 = sns.jointplot(data=data, x=target_feature, y=i, hue='sample')
        plt.savefig('%scomparison_of_%s_vs_%s.pdf' % (master_folder, target_feature, i))
        plt.close()"""

# radial distribution
"""feature = 'radial_distribution_relative_r'
data_DM[feature] = [dat.str_to_float(data_DM[feature][i]) for i in range(len(data_DM))]
data_HSR[feature] = [dat.str_to_float(data_HSR[feature][i]) for i in range(len(data_HSR))]

number_nuclear_DM = len(data_DM)
number_nuclear_HSR = len(data_HSR)

FISH_mean_curve_DM, FISH_ci_lower_DM, FISH_ci_higher_DM = dat.mean_list(data_DM[feature].tolist())
FISH_mean_curve_HSR, FISH_ci_lower_HSR, FISH_ci_higher_HSR = dat.mean_list(data_HSR[feature].tolist())

r = np.arange(0, 1, relative_radial_interval)

plt.subplots(figsize=(6, 4))
for i in range(len(data_HSR)):
    plt.plot(r, data_HSR[feature][i], alpha=0.05, color=[colors[0][j]+0.1 for j in range(len(colors[0]))])
for i in range(len(data_DM)):
    plt.plot(r, data_DM[feature][i], alpha=0.05, color=[colors[1][j]+0.1 for j in range(len(colors[1]))])
plt.plot(r, FISH_mean_curve_DM, color=colors[1], label='DM, n=%s' % number_nuclear_DM)
plt.plot(r, FISH_ci_lower_DM, color=colors[1], linestyle='--', linewidth=0.5)
plt.plot(r, FISH_ci_higher_DM, color=colors[1], linestyle='--', linewidth=0.5)
plt.plot(r, FISH_mean_curve_HSR, color=colors[0], label='HSR, n=%s' % number_nuclear_HSR)
plt.plot(r, FISH_ci_lower_HSR, color=colors[0], linestyle='--', linewidth=0.5)
plt.plot(r, FISH_ci_higher_HSR, color=colors[0], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='#FF4500', linestyle='--')
plt.xlabel('relative r')
plt.ylabel('normalized distribution')
plt.ylim([0, 2])
plt.legend()
plt.savefig('%s/%s_comparison.pdf' % (master_folder, feature))
plt.close()

plt.subplots(figsize=(6, 4))
for i in range(len(data_HSR)):
    plt.plot(r, data_HSR[feature][i], alpha=0.05, color=[colors[0][j]+0.1 for j in range(len(colors[0]))])
plt.plot(r, FISH_mean_curve_HSR, color=colors[0], label='HSR, n=%s' % number_nuclear_HSR)
plt.plot(r, FISH_ci_lower_HSR, color=colors[0], linestyle='--', linewidth=0.5)
plt.plot(r, FISH_ci_higher_HSR, color=colors[0], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='#FF4500', linestyle='--')
plt.xlabel('relative r')
plt.ylabel('normalized distribution')
plt.ylim([0, 2])
plt.legend()
plt.savefig('%s/%s_HSR.pdf' % (master_folder, feature))
plt.close()

plt.subplots(figsize=(6, 4))
for i in range(len(data_DM)):
    plt.plot(r, data_DM[feature][i], alpha=0.05, color=[colors[1][j]+0.1 for j in range(len(colors[1]))])
plt.plot(r, FISH_mean_curve_DM, color=colors[1], label='DM, n=%s' % number_nuclear_DM)
plt.plot(r, FISH_ci_lower_DM, color=colors[1], linestyle='--', linewidth=0.5)
plt.plot(r, FISH_ci_higher_DM, color=colors[1], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='#FF4500', linestyle='--')
plt.xlabel('relative r')
plt.ylabel('normalized distribution')
plt.ylim([0, 2])
plt.legend()
plt.savefig('%s/%s_DM.pdf' % (master_folder, feature))
plt.close()"""

# radial distribution vs intensity/g_value
"""data_DM['sample'] = ['DM'] * len(data_DM)
data_HSR['sample'] = ['HSR'] * len(data_HSR)

feature = 'radial_distribution_relative_r'
data_DM[feature] = [dat.str_to_float(data_DM[feature][i]) for i in range(len(data_DM))]
data_HSR[feature] = [dat.str_to_float(data_HSR[feature][i]) for i in range(len(data_HSR))]

data_DM['g'] = [dat.str_to_float(data_DM['g'][i]) for i in range(len(data_DM))]
data_HSR['g'] = [dat.str_to_float(data_HSR['g'][i]) for i in range(len(data_HSR))]

data = pd.concat([data_HSR, data_DM], axis=0, ignore_index=True)
data['g_value'] = [(data['g'][i][1] + data['g'][i][2] + data['g'][i][3] + data['g'][i][4] + data['g'][i][5])*0.2
                   for i in range(len(data))]

data['R4'] = [(data[feature][i][2] + data[feature][i][3] + data[feature][i][4] + data[feature][i][5] +
               data[feature][i][6])*0.2 for i in range(len(data))]
data['R66'] = [(data[feature][i][64] + data[feature][i][65] + data[feature][i][66] + data[feature][i][67] +
                data[feature][i][68])*0.2 for i in range(len(data))]
data['R4/R66'] = data['R4']/data['R66']

data['MYC_intensity/FISH_intensity'] = data['MYC_total_intensity']/data['FISH_total_intensity_nuclear']

sns.set_palette(sns.color_palette(colors))

# ax1 = sns.jointplot(data=data, x='R35l', y='MYC_intensity/FISH_intensity', hue='sample')
ax1 = sns.jointplot(data=data, x='R35l', y='g_value', hue='sample')

_, g_fit_r2_DM, g_fit_a_DM, g_fit_b_DM = \
    mat.fitting_linear(data[data['sample'] == 'DM']['R35l'].tolist(),
                          data[data['sample'] == 'DM']['MYC_intensity/FISH_intensity'].tolist())
_, g_fit_r2_HSR, g_fit_a_HSR, g_fit_b_HSR = \
    mat.fitting_linear(data[data['sample'] == 'HSR']['R35l'].tolist(),
                          data[data['sample'] == 'HSR']['MYC_intensity/FISH_intensity'].tolist())

x = np.arange(0, 0.2, 0.01)
y_DM = g_fit_a_DM * x + g_fit_b_DM
y_HSR = g_fit_a_HSR * x + g_fit_b_HSR
ax1.ax_joint.plot(x, y_HSR, linewidth=2, color=colors[0], linestyle='--')
ax1.ax_joint.plot(x, y_DM, linewidth=2, color=colors[1], linestyle='--')

# plt.savefig('%scomparison_of_radial_distribution_vs_intensity_residue.pdf' % master_folder)
plt.savefig('%scomparison_of_radial_distribution_vs_g_value.pdf' % master_folder)
plt.close()"""

# direction map of ecDNA centroid
"""data_DM['sample'] = ['DM'] * len(data_DM)
data_HSR['sample'] = ['HSR'] * len(data_HSR)

data_DM['ecDNA_localization_from_centroid'] = \
    [dat.str_to_list_of_float(data_DM['ecDNA_localization_from_centroid'][i], 2) for i in range(len(data_DM))]
data_HSR['ecDNA_localization_from_centroid'] = \
    [dat.str_to_list_of_float(data_HSR['ecDNA_localization_from_centroid'][i], 2) for i in range(len(data_HSR))]

x_DM, y_DM = dat.list_separation(data_DM, 'ecDNA_localization_from_centroid')
x_HSR, y_HSR = dat.list_separation(data_HSR, 'ecDNA_localization_from_centroid')

temp = pd.DataFrame({'x': x_HSR+x_DM, 'y': y_HSR+y_DM, 'sample': ['HSR']*len(x_HSR)+['DM']*len(x_DM)})

colors = [(0.85, 0.35, 0.25), (0.2, 0.2, 0.2)]
sns.set_palette(sns.color_palette(colors))
ax1 = sns.displot(data=temp[temp['sample'] == 'DM'], x='x', y='y', bins=15)
plt.ylim([-100, 100])
plt.xlim([-100, 100])
plt.savefig('%s/direction_map_of_ecDNA_centroid_DM.pdf' % master_folder)
plt.close()
colors = [(0.2, 0.2, 0.2), (0.85, 0.35, 0.25)]
sns.set_palette(sns.color_palette(colors))
ax2 = sns.displot(data=temp[temp['sample'] == 'HSR'], x='x', y='y', bins=15)
plt.ylim([-100, 100])
plt.xlim([-100, 100])
plt.savefig('%s/direction_map_of_ecDNA_centroid_HSR.pdf' % master_folder)
plt.close()"""

# auto-correlation, different int_thresh, different seed number
data_auto['g'] = [dat.str_to_float(data_auto['g'][i]) for i in range(len(data_auto))]
data_auto['nan'] = [data_auto['g'][i][0] for i in range(len(data_auto))]
data_auto.dropna(subset=['nan'], inplace=True)
data_auto.reset_index(drop=True, inplace=True)

df = data_auto.pivot(columns='int_thresh', index='dots', values='g_value')
# ax1 = sns.heatmap(data=df, annot=True, cmap='viridis', cbar=True, cbar_kws={'label': 'g_value'}, square=True)
ax1 = sns.heatmap(data=df, cmap='viridis', cbar=True, cbar_kws={'label': 'g_value'}, square=True)
plt.savefig('%s/effect_of_dots_and_intThresh_on_auto_correlation_fov0_i2.pdf' % master_folder)
plt.close()

