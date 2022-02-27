import pandas as pd
import numpy as np

# input parameters
exp_name = '20220204_CtrlAndJQ1'
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220204_CtrlAndJQ1/"
prefix = 'ctrl'
interval = 0.01
total_FOV = 20

data_curve = pd.DataFrame(columns=['FOV', 'nuclear', 'interval', 'density_FISH', 'density_nuclear'])

for fov in range(total_FOV):
    print("Start analyzing FOV %s/%s..." % (fov+1, total_FOV))
    data_r = pd.read_csv('%sdistance_analysis/%s/summary_r_%s_FOV%s.txt' % (master_folder, prefix, prefix, fov), na_values=['.'], sep='\t')

    for i in range(int(max(data_r['nuclear']))):
        data_temp = data_r[data_r['nuclear'] == i + 1]
        if len(data_temp) != 0:
            data_temp = data_temp.sort_values(by='relative_r')
            y_FISH = []
            y_nucleus = []
            for j in np.arange(0, 1, interval):
                if j != 1 - interval:
                    data_temp2 = data_temp[(data_temp['relative_r'] >= j) & (data_temp['relative_r'] < j + interval)]
                else:
                    data_temp2 = data_temp[data_temp['relative_r'] >= j]
                if (len(data_temp2)) != 0:
                    y_FISH.append(sum(data_temp2['intensity_FISH'].tolist()) / sum(data_temp2['mean_intensity_FISH'].tolist()))
                    y_nucleus.append(sum(data_temp2['intensity_nucleus'].tolist()) / sum(data_temp2['mean_intensity_nucleus'].tolist()))
                else:
                    y_FISH.append(0)
                    y_nucleus.append(0)
            data_curve.loc[len(data_curve.index)] = [fov, i + 1, interval, y_FISH, y_nucleus]

data_curve.to_csv('%ssummary_curve_%s.txt' % (master_folder, prefix), index=False, sep='\t')

print("END!")
