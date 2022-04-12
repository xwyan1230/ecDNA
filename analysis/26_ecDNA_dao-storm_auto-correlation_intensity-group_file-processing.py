import pandas as pd
import numpy as np
import shared.dataframe as dat
import matplotlib.pyplot as plt

# input parameters
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220301_ecDNA_ctrlAndJQ1_NatashaFile/"

data = pd.read_csv('%sauto_correlation_JQ13hr_1.txt' % master_folder, na_values=['.'], sep='\t')
data['g_FISH'] = [dat.str_to_float(data['g_FISH'][i]) for i in range(len(data))]
data['nan'] = [data['g_FISH'][i][0] for i in range(len(data))]
data.dropna(subset=['nan'], inplace=True)
data.reset_index(drop=True, inplace=True)
data['g_nuclear'] = [dat.str_to_float(data['g_nuclear'][i]) for i in range(len(data))]
data['nuclear_centroid'] = [dat.str_to_float(data['nuclear_centroid'][i]) for i in range(len(data))]
data['FISH_int'] = data['nuclear_area']*data['FISH_mean_int']
data['centroid_x'] = [data['nuclear_centroid'][i][0] for i in range(len(data))]
data['centroid_y'] = [data['nuclear_centroid'][i][1] for i in range(len(data))]

z = int(max(data['z'])/2)
data_z = data[data['z'] == z]

data_nuclear = pd.DataFrame(columns=['nuclear', 'z', 'z_min', 'z_max', 'nuclear_area', 'FISH_int', 'g_FISH', 'g_nuclear'])
data_nuclear_mean = pd.DataFrame(columns=['nuclear', 'FISH_mean_int', 'g_FISH', 'g_nuclear'])

for i in range(len(data_z)):
    data_temp = data[(data['centroid_x'] > data_z['centroid_x'].tolist()[i]-50) &
                     (data['centroid_x'] < data_z['centroid_x'].tolist()[i]+50) &
                     (data['centroid_y'] > data_z['centroid_y'].tolist()[i]-50) &
                     (data['centroid_y'] < data_z['centroid_y'].tolist()[i]+50)]
    if len(data_temp) == len(set(data_temp['z'].tolist())):
        data_temp = data_temp.sort_values(by='FISH_int', ascending=False)
        data_temp.reset_index(drop=True, inplace=True)
        z_min = min(data_temp['z'].tolist())
        z_max = max(data_temp['z'].tolist())
        data_nuclear.loc[len(data_nuclear.index)] = [i, data_temp['z'].tolist()[0], z_min, z_max,
                                                     data_temp['nuclear_area'].tolist()[0],
                                                     data_temp['FISH_int'].tolist()[0], data_temp['g_FISH'].tolist()[0],
                                                     data_temp['g_nuclear'].tolist()[0]]
        data_nuclear.loc[len(data_nuclear.index)] = [i, data_temp['z'].tolist()[1], z_min, z_max,
                                                     data_temp['nuclear_area'].tolist()[1],
                                                     data_temp['FISH_int'].tolist()[1], data_temp['g_FISH'].tolist()[1],
                                                     data_temp['g_nuclear'].tolist()[1]]
        data_nuclear.loc[len(data_nuclear.index)] = [i, data_temp['z'].tolist()[2], z_min, z_max,
                                                     data_temp['nuclear_area'].tolist()[2],
                                                     data_temp['FISH_int'].tolist()[2], data_temp['g_FISH'].tolist()[2],
                                                     data_temp['g_nuclear'].tolist()[2]]
        g_FISH_mean = list((np.array(data_temp['g_FISH'].tolist()[0]) + np.array(data_temp['g_FISH'].tolist()[1]) +
                            np.array(data_temp['g_FISH'].tolist()[2]))*1.0/3)
        g_nuclear_mean = list((np.array(data_temp['g_nuclear'].tolist()[0]) + np.array(data_temp['g_nuclear'].tolist()[1])
                               + np.array(data_temp['g_nuclear'].tolist()[2])) * 1.0 / 3)
        FISH_int_mean = (data_temp['FISH_mean_int'][0] + data_temp['FISH_mean_int'][1] + data_temp['FISH_mean_int'][2]) * 1.0 / 3
        data_nuclear_mean.loc[len(data_nuclear_mean.index)] = [i, FISH_int_mean, g_FISH_mean, g_nuclear_mean]
    else:
        print("Nuclei too close.")

data_nuclear.to_csv('%sauto_correlation_nuclear_JQ13hr_1.txt' % master_folder, index=False, sep='\t')
data_nuclear_mean.to_csv('%sauto_correlation_nuclear_mean_JQ13hr_1.txt' % master_folder, index=False, sep='\t')

print("END!")