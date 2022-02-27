import pandas as pd
import numpy as np
import shared.dataframe as dat
import matplotlib.pyplot as plt

# input parameters
exp_name = '20220204_CtrlAndJQ1'
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220204_CtrlAndJQ1/"
prefix = ['ctrl', 'JQ1']
color = ['#FFA500', '#40E0D0']
rmax = 60
mean_curve = []
ci_lower = []
ci_higher = []
n_nuclear = []

for p in prefix:
    data = pd.read_csv('%sauto_correlation_%s.txt' % (master_folder, p), na_values=['.'], sep='\t')
    data['g_FISH'] = [dat.str_to_float(data['g_FISH'][i]) for i in range(len(data))]
    data['dg_FISH'] = [dat.str_to_float(data['dg_FISH'][i]) for i in range(len(data))]
    data['g_nuclear'] = [dat.str_to_float(data['g_nuclear'][i]) for i in range(len(data))]
    data['dg_nuclear'] = [dat.str_to_float(data['dg_nuclear'][i]) for i in range(len(data))]
    data['nan'] = [data['g_FISH'][i][0] for i in range(len(data))]
    data.dropna(subset=['nan'], inplace=True)
    data.reset_index(drop=True, inplace=True)

    number_nuclear = len(data)
    n_nuclear.append(number_nuclear)

    nuclear_mean_curve, _, _ = dat.mean_list(data['g_nuclear'].tolist())

    data['g_FISH_adjust'] = [np.array(data['g_FISH'][i]) / np.array(nuclear_mean_curve) for i in range(len(data))]
    FISH_mean_curve, FISH_ci_lower, FISH_ci_higher = dat.mean_list(data['g_FISH_adjust'].tolist())
    mean_curve.append(FISH_mean_curve)
    ci_lower.append(FISH_ci_lower)
    ci_higher.append(FISH_ci_higher)
    
    r = np.arange(0, rmax+1, 1)

    plt.subplots(figsize=(6, 4))
    for i in range(len(data)):
        plt.plot(r, data['g_FISH_adjust'][i], alpha=0.1, color='#FFD700')
    plt.plot(r, FISH_mean_curve, color='#FF4500', label='mean, n=%s' % number_nuclear)
    plt.plot(r, FISH_ci_lower, color='#FF4500', linestyle='--', linewidth=0.5)
    plt.plot(r, FISH_ci_higher, color='#FF4500', linestyle='--', linewidth=0.5)
    plt.xlabel('r')
    plt.ylabel('g')
    plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
    plt.savefig('%s/auto_correlation_%s.pdf' % (master_folder, p))
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
plt.savefig('%s/auto_correlation_comparison.pdf' % master_folder)
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
    plt.savefig('%s/auto_correlation_ci_%s.pdf' % (master_folder, prefix[i]))
    plt.close()


