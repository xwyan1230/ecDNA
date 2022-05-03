import pandas as pd
import numpy as np
import shared.dataframe as dat
import matplotlib.pyplot as plt

# input parameters
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/211102_3hr_JQ1washout/"
prefix = ['DMSO3hr', 'JQ13hr', '3hrJQ1_3hrDMSO_WO', '3hrJQ1_3hr1uMtriptolide_WO']
color = ['#FFA500', '#40E0D0', '#DA70D6', '#00CED1']
pixel_size = 40  # nm (Zeiss confocal scope)
cell_avg_size = 10  # um (Colo)
local_size = int(0.9 * cell_avg_size * 1000/pixel_size)  # single FOV only, assuming squared pixel, ~150
rmax = int(0.67 * local_size)  # ~100
mean_curve = []
ci_lower = []
ci_higher = []
n_nuclear = []
feature = 'g_value'

for p in prefix:
    data = pd.read_csv('%s%s_different_int_threshold_g.txt' % (master_folder, p), na_values=['.'], sep='\t')
    data[feature] = [dat.str_to_float(data[feature][i]) for i in range(len(data))]

    number_nuclear = len(data)
    n_nuclear.append(number_nuclear)

    FISH_mean_curve, FISH_ci_lower, FISH_ci_higher = dat.mean_list(data[feature].tolist())
    mean_curve.append(FISH_mean_curve)
    ci_lower.append(FISH_ci_lower)
    ci_higher.append(FISH_ci_higher)

    # r = np.arange(0, rmax + 1, 1)
    r = np.arange(0, 2000, 100)

    plt.subplots(figsize=(6, 4))
    for i in range(len(data)):
        plt.plot(r, data[feature][i], alpha=0.2, color='#FFD700')
    plt.plot(r, FISH_mean_curve, color='#FF4500', label='mean, n=%s' % number_nuclear)
    plt.plot(r, FISH_ci_lower, color='#FF4500', linestyle='--', linewidth=0.5)
    plt.plot(r, FISH_ci_higher, color='#FF4500', linestyle='--', linewidth=0.5)
    plt.xlabel('r')
    plt.ylabel('g')
    plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
    plt.savefig('%s/g_value_%s.pdf' % (master_folder, p))
    plt.close()

plt.subplots(figsize=(6, 4))
for i in range(len(prefix)):
    plt.plot(r, mean_curve[i], color=color[i], label='%s, n=%s' % (prefix[i], n_nuclear[i]))
    plt.plot(r, ci_lower[i], color=color[i], linestyle='--', linewidth=0.5)
    plt.plot(r, ci_higher[i], color=color[i], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='#FF4500', linestyle='--')
plt.xlabel('r')
plt.ylabel('g')
plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
plt.savefig('%s/g_value_comparison.pdf' % master_folder)
plt.close()

for i in range(len(prefix)):
    plt.subplots(figsize=(6, 4))
    plt.plot(r, mean_curve[i], color=color[i], label='%s, n=%s' % (prefix[i], n_nuclear[i]))
    plt.plot(r, ci_lower[i], color=color[i], linestyle='--', linewidth=0.5)
    plt.plot(r, ci_higher[i], color=color[i], linestyle='--', linewidth=0.5)
    plt.axhline(y=1, color='#FF4500', linestyle='--')
    plt.xlabel('r')
    plt.ylabel('g')
    plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
    plt.savefig('%s/g_value_ci_%s.pdf' % (master_folder, prefix[i]))
    plt.close()