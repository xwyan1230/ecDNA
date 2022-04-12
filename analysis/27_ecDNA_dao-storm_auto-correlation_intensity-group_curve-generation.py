import pandas as pd
import numpy as np
import shared.dataframe as dat
import matplotlib.pyplot as plt

# input parameters
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220301_ecDNA_ctrlAndJQ1_NatashaFile/"
color = ['#FFA500', '#40E0D0', '#DA70D6', '#00CED1']
rmax = 100

data = pd.read_csv('%sauto_correlation_nuclear_mean_JQ13hr_1.txt' % master_folder, na_values=['.'], sep='\t')
data['g_FISH'] = [dat.str_to_float(data['g_FISH'][i]) for i in range(len(data))]
data['g_nuclear'] = [dat.str_to_float(data['g_nuclear'][i]) for i in range(len(data))]

number_nuclear = len(data)

FISH_mean_curve, _, _ = dat.mean_list(data['g_FISH'].tolist())
nuclear_mean_curve, _, _ = dat.mean_list(data['g_nuclear'].tolist())

r = np.arange(0, rmax + 1, 1)

plt.subplots(figsize=(6, 4))
for i in range(len(data)):
    plt.plot(r, data['g_FISH'][i], alpha=0.2, color='#FFD700')
    plt.plot(r, data['g_nuclear'][i], alpha=0.2, color='#40E0D0')
plt.plot(r, FISH_mean_curve, color='#FF4500', label='FISH mean, n=%s' % number_nuclear)
plt.plot(r, nuclear_mean_curve, color='#00CED1', label='nuclear mean, n=%s' % number_nuclear)
plt.xlabel('r')
plt.ylabel('g')
plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
plt.savefig('%s/auto_correlation_JQ13hr_1.pdf' % master_folder)
plt.ylim([-0.5, 20.5])
plt.savefig('%s/auto_correlation_JQ13hr_1_part.pdf' % master_folder)
plt.close()

"""data = pd.read_csv('%sauto_correlation_nuclear_mean_JQ13hr_1.txt' % master_folder, na_values=['.'], sep='\t')

data['FISH_tip'] = [data['g_FISH'][i][5] for i in range(len(data))]
data['nuclear_tip'] = [data['g_nuclear'][i][5] for i in range(len(data))]

plt.subplots(figsize=(6, 4))
plt.scatter(data['FISH_mean_int'], data['FISH_tip'])
plt.xlabel('FISH_int_mean')
plt.ylabel('auto_correlation_peak')
plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
plt.savefig('%s/auto_correlation_intensity-scatter_JQ13hr_1.pdf' % master_folder)
plt.close()"""