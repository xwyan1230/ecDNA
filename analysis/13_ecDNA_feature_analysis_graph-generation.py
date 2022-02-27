import pandas as pd
import numpy as np
import shared.dataframe as dat
import matplotlib.pyplot as plt
import seaborn as sns

# input parameters
exp_name = '20220204_CtrlAndJQ1'
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220204_CtrlAndJQ1/"
prefix = ['ctrl', 'JQ1']
color = ['#FFA500', '#40E0D0']
ecDNA_number = []
ecDNA_area = []
ecDNA_pc = []
features = ['ecDNA_number', 'ecDNA_area', 'ecDNA_pc']
bin = [np.arange(0, 55, 0.5), np.arange(0, 500, 5), np.arange(0, 1, 0.01)]

for p in prefix:
    data = pd.read_csv('%sfeature_analysis_%s.txt' % (master_folder, p), na_values=['.'], sep='\t')
    ecDNA_number_temp = data['ecDNA_number'].tolist()
    ecDNA_number.append(ecDNA_number_temp)

    data = data.loc[data['ecDNA_number'] != 0]
    data.reset_index(drop=True, inplace=True)
    data['ecDNA_area'] = [dat.str_to_float(data['ecDNA_area'][i]) for i in range(len(data))]
    data['ecDNA_mean_int'] = [dat.str_to_float(data['ecDNA_mean_int'][i]) for i in range(len(data))]

    ecDNA_area_temp = []
    for i in range(len(data)):
        ecDNA_area_temp = ecDNA_area_temp + data['ecDNA_area'][i]
    ecDNA_area.append(ecDNA_area_temp)

    data['ecDNA_int'] = [np.array(data['ecDNA_area'][i]) * np.array(data['ecDNA_mean_int'][i]) for i in range(len(data))]
    data['ecDNA_int_sum'] = [sum(data['ecDNA_int'][i]) for i in range(len(data))]
    data['pc'] = data['ecDNA_int_sum']*1.0/(data['nuclear_area']*data['nuclear_mean_int'])
    ecDNA_pc_temp = data['pc'].tolist()
    ecDNA_pc.append(ecDNA_pc_temp)

df = pd.DataFrame({'sample': prefix, 'ecDNA_number': ecDNA_number, 'ecDNA_area': ecDNA_area, 'ecDNA_pc': ecDNA_pc})

for j in range(len(features)):
    plt.subplots(figsize=(6, 4))
    for i in range(len(prefix)):
        n = len(df[features[j]][i])
        weights = np.ones_like(df[features[j]][i]) / len(df[features[j]][i])
        plt.hist(df[features[j]][i], weights=weights, bins=bin[j], color=color[i], edgecolor='w', alpha=0.5, label='%s, n=%s' % (prefix[i], n))
    plt.xlabel('%s' % features[j])
    plt.ylabel('AU')
    plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
    plt.savefig('%s/%s_comparison.pdf' % (master_folder, features[j]))
    plt.close()
