import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# input parameters
exp_name = '20220204_CtrlAndJQ1'
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220204_CtrlAndJQ1/"
interval = 0.01

data_r = pd.read_csv('%ssummary_r.txt' % master_folder, na_values=['.'], sep='\t')

for i in range(int(max(data_r['nuclear']))):
    data_temp = data_r[data_r['nuclear'] == i + 1]
    data_temp = data_temp.sort_values(by='relative_r')
    x = []
    y_FISH = []
    y_nucleus = []
    for j in np.arange(0, 1, interval):
        if j != 1 - interval:
            data_temp2 = data_temp[(data_temp['relative_r'] >= j) & (data_temp['relative_r'] < j + interval)]
        else:
            data_temp2 = data_temp[data_temp['relative_r'] >= j]
        if (len(data_temp2)) != 0:
            x.append(j + 0.5 * interval)
            y_FISH.append(sum(data_temp2['intensity_FISH'].tolist()) / sum(data_temp2['mean_intensity_FISH'].tolist()))
            y_nucleus.append(sum(data_temp2['intensity_nucleus'].tolist()) / sum(data_temp2['mean_intensity_nucleus'].tolist()))

    plt.subplots(figsize=(6, 4))
    plt.plot(x, y_FISH, color=(0.85, 0.35, 0.25), label='nuclear %s ecDNA' % (i + 1))
    plt.plot(x, y_nucleus, color='#1E90FF', label='nuclear %s chromosome' % (i + 1))
    plt.xlabel('relative r')
    plt.ylabel('intensity density')
    plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
    plt.savefig('%s/nuclear%s.pdf' % (master_folder, i + 1))
    plt.close()

