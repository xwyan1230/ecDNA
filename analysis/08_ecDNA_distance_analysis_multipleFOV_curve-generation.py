import pandas as pd
import numpy as np
import shared.dataframe as dat
import matplotlib.pyplot as plt

# input parameters
exp_name = '20220204_CtrlAndJQ1'
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220204_CtrlAndJQ1/"
prefix = ['ctrl', 'JQ1']
mean_curve = []
n_nuclear = []

for p in prefix:
    data_curve = pd.read_csv('%ssummary_curve_%s.txt' % (master_folder, p), na_values=['.'], sep='\t')
    data_curve['density_FISH'] = [dat.str_to_float(data_curve['density_FISH'][i]) for i in range(len(data_curve))]
    data_curve['density_nuclear'] = [dat.str_to_float(data_curve['density_nuclear'][i]) for i in range(len(data_curve))]
    interval = data_curve['interval'][0]
    number_nuclear = len(data_curve)
    n_nuclear.append(number_nuclear)

    nuclear_mean_curve = dat.mean_list(data_curve['density_nuclear'].tolist())
    FISH_mean_curve = dat.mean_list(data_curve['density_FISH'].tolist())
    FISH_adjust_curve = np.array(FISH_mean_curve)/np.array(nuclear_mean_curve)
    mean_curve.append(FISH_adjust_curve)

    data_curve['density_FISH_adjust'] = [np.array(data_curve['density_FISH'][i])/np.array(nuclear_mean_curve) for i in range(len(data_curve))]

    x = np.arange(0, 1, interval)

    plt.subplots(figsize=(6, 4))
    for i in range(len(data_curve)):
        x_exclude, y_exclude = dat.list_exclude_zero(x, data_curve['density_FISH_adjust'][i])
        plt.plot(x_exclude, y_exclude, alpha=0.1, color='#1E90FF')
    plt.plot(x, FISH_adjust_curve, color=(0.85, 0.35, 0.25), label='mean, n=%s' % number_nuclear)
    plt.xlabel('Relative distance from centroid (AU)')
    plt.ylabel('Normalized intensity density (AU)')
    plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
    plt.savefig('%s/curve_%s.pdf' % (master_folder, p))
    plt.close()

plt.subplots(figsize=(6, 4))
for i in range(len(prefix)):
    plt.plot(x, mean_curve[i], label='%s, n=%s' % (prefix[i], n_nuclear[i]))
plt.axhline(y=1, color='r', linestyle='--')
plt.xlabel('Relative distance from centroid (AU)')
plt.ylabel('Normalized intensity density (AU)')
plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
plt.savefig('%s/curve_comparison.pdf' % master_folder)
plt.close()



