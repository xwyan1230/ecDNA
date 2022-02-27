import pandas as pd
import numpy as np
import shared.dataframe as dat
import matplotlib.pyplot as plt

# input parameters
exp_name = '20220204_CtrlAndJQ1'
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220204_CtrlAndJQ1/"
prefix = ['ctrl', 'JQ1']
color = ['#FFA500', '#40E0D0']
mean_curve = []
ci_lower = []
ci_higher = []
n_nuclear = []

for p in prefix:
    data_curve = pd.read_csv('%ssummary_curve_%s.txt' % (master_folder, p), na_values=['.'], sep='\t')
    data_curve['density_FISH'] = [dat.str_to_float(data_curve['density_FISH'][i]) for i in range(len(data_curve))]
    data_curve['density_nuclear'] = [dat.str_to_float(data_curve['density_nuclear'][i]) for i in range(len(data_curve))]
    interval = data_curve['interval'][0]
    number_nuclear = len(data_curve)
    n_nuclear.append(number_nuclear)

    nuclear_mean_curve, _, _ = dat.mean_list(data_curve['density_nuclear'].tolist())
    data_curve['density_FISH_adjust'] = [np.array(data_curve['density_FISH'][i]) / np.array(nuclear_mean_curve) for i in
                                         range(len(data_curve))]
    FISH_mean_curve, FISH_ci_lower, FISH_ci_higher = dat.mean_list(data_curve['density_FISH_adjust'].tolist())
    mean_curve.append(FISH_mean_curve)
    ci_lower.append(FISH_ci_lower)
    ci_higher.append(FISH_ci_higher)

    x = np.arange(0, 1, interval)

    plt.subplots(figsize=(6, 4))
    for i in range(len(data_curve)):
        x_exclude, y_exclude = dat.list_exclude_zero(x, data_curve['density_FISH_adjust'][i])
        plt.plot(x_exclude, y_exclude, alpha=0.1, color='#FFD700')
    x_exclude, y_exclude = dat.list_exclude_zero(x, FISH_mean_curve)
    plt.plot(x_exclude, y_exclude, color='#FF4500', label='mean, n=%s' % number_nuclear)
    x_exclude, y_exclude = dat.list_exclude_zero(x, FISH_ci_lower)
    plt.plot(x_exclude, y_exclude, color='#FF4500', linestyle='--', linewidth=0.5)
    x_exclude, y_exclude = dat.list_exclude_zero(x, FISH_ci_higher)
    plt.plot(x_exclude, y_exclude, color='#FF4500', linestyle='--', linewidth=0.5)
    plt.xlabel('Relative distance from centroid (AU)')
    plt.ylabel('Normalized intensity density (AU)')
    plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
    plt.savefig('%s/curve_%s.pdf' % (master_folder, p))
    plt.close()

plt.subplots(figsize=(6, 4))
for i in range(len(prefix)):
    x_exclude, y_exclude = dat.list_exclude_zero(x, mean_curve[i])
    plt.plot(x_exclude, y_exclude, color = color[i], label='%s, n=%s' % (prefix[i], n_nuclear[i]))
    x_exclude, y_exclude = dat.list_exclude_zero(x, ci_lower[i])
    plt.plot(x_exclude, y_exclude, color = color[i], linestyle='--', linewidth=0.5)
    x_exclude, y_exclude = dat.list_exclude_zero(x, ci_higher[i])
    plt.plot(x_exclude, y_exclude, color = color[i], linestyle='--', linewidth=0.5)
plt.axhline(y=1, color='#FF4500', linestyle='--')
plt.xlabel('Relative distance from centroid (AU)')
plt.ylabel('Normalized intensity density (AU)')
plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
plt.savefig('%s/curve_comparison.pdf' % master_folder)
plt.close()

for i in range(len(prefix)):
    plt.subplots(figsize=(6, 4))
    x_exclude, y_exclude = dat.list_exclude_zero(x, mean_curve[i])
    plt.plot(x_exclude, y_exclude, color = color[i], label='%s, n=%s' % (prefix[i], n_nuclear[i]))
    x_exclude, y_exclude = dat.list_exclude_zero(x, ci_lower[i])
    plt.plot(x_exclude, y_exclude, color = color[i], linestyle='--', linewidth=0.5)
    x_exclude, y_exclude = dat.list_exclude_zero(x, ci_higher[i])
    plt.plot(x_exclude, y_exclude, color = color[i], linestyle='--', linewidth=0.5)
    plt.axhline(y=1, color='#FF4500', linestyle='--')
    plt.xlabel('Relative distance from centroid (AU)')
    plt.ylabel('Normalized intensity density (AU)')
    plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
    plt.savefig('%s/curve_ci_%s.pdf' % (master_folder, prefix[i]))
    plt.close()



