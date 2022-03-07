import pandas as pd
import numpy as np
import shared.dataframe as dat

# input parameters
exp_name = '20220204_CtrlAndJQ1'
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220204_CtrlAndJQ1_zstack/"
sample = 'ctrl'
total_FOV = 20

data = pd.read_csv('%sauto_correlation_%s.txt' % (master_folder, sample), na_values=['.'], sep='\t')
data['g_FISH'] = [dat.str_to_float(data['g_FISH'][i]) for i in range(len(data))]
data['nan'] = [data['g_FISH'][i][0] for i in range(len(data))]
data.dropna(subset=['nan'], inplace=True)
data.reset_index(drop=True, inplace=True)
data['nuclear_centroid'] = [dat.str_to_float(data['nuclear_centroid'][i]) for i in range(len(data))]
data['FISH_int'] = data['nuclear_area']*data['FISH_mean_int']
data['centroid_x'] = [data['nuclear_centroid'][i][0] for i in range(len(data))]
data['centroid_y'] = [data['nuclear_centroid'][i][1] for i in range(len(data))]

data_nuclear = pd.DataFrame(columns=['FOV', 'nuclear', 'z', 'nuclear_area', 'FISH_int', 'centroid_x', 'centroid_y', 'g_FISH'])
data_nuclear_mean = pd.DataFrame(columns=['FOV', 'nuclear', 'centroid_x', 'centroid_y', 'g_FISH'])

for fov in range(total_FOV):
    data_fov = data[data['FOV'] == fov]
    z = int(max(data_fov['z'])/2)
    data_z = data_fov[data_fov['z'] == z]

    for i in range(len(data_z)):
        data_temp = data[(data['centroid_x'] > data_z['centroid_x'].tolist()[i]-10) &
                         (data['centroid_x'] < data_z['centroid_x'].tolist()[i]+10) &
                         (data['centroid_y'] > data_z['centroid_y'].tolist()[i]-10) &
                         (data['centroid_y'] < data_z['centroid_y'].tolist()[i]+10)]
        if (len(data_temp) == len(set(data_temp['z'].tolist()))) & (len(data_temp) > 2):
            data_temp = data_temp.sort_values(by='FISH_int', ascending=False)
            data_temp.reset_index(drop=True, inplace=True)
            data_nuclear.loc[len(data_nuclear.index)] = [fov, i + 1, data_temp['z'].tolist()[0],
                                                         data_temp['nuclear_area'].tolist()[0],
                                                         data_temp['FISH_int'].tolist()[0],
                                                         data_temp['centroid_x'].tolist()[0],
                                                         data_temp['centroid_y'].tolist()[0],
                                                         data_temp['g_FISH'].tolist()[0]]
            data_nuclear.loc[len(data_nuclear.index)] = [fov, i + 1, data_temp['z'].tolist()[1],
                                                         data_temp['nuclear_area'].tolist()[1],
                                                         data_temp['FISH_int'].tolist()[1],
                                                         data_temp['centroid_x'].tolist()[1],
                                                         data_temp['centroid_y'].tolist()[1],
                                                         data_temp['g_FISH'].tolist()[1]]
            data_nuclear.loc[len(data_nuclear.index)] = [fov, i + 1, data_temp['z'].tolist()[2],
                                                         data_temp['nuclear_area'].tolist()[2],
                                                         data_temp['FISH_int'].tolist()[2],
                                                         data_temp['centroid_x'].tolist()[2],
                                                         data_temp['centroid_y'].tolist()[2],
                                                         data_temp['g_FISH'].tolist()[2]]
            g_FISH_mean = list((np.array(data_temp['g_FISH'].tolist()[0]) + np.array(data_temp['g_FISH'].tolist()[1]) +
                           np.array(data_temp['g_FISH'].tolist()[2]))*1.0/3)
            data_nuclear_mean.loc[len(data_nuclear_mean.index)] = [fov, i+1, data_temp['centroid_x'].tolist()[0],
                                                                   data_temp['centroid_y'].tolist()[0], g_FISH_mean]
        else:
            print("Nuclei too close.")

data_nuclear.to_csv('%sauto_correlation_nuclear_%s.txt' % (master_folder, sample), index=False, sep='\t')
data_nuclear_mean.to_csv('%sauto_correlation_nuclear_mean_%s.txt' % (master_folder, sample), index=False, sep='\t')

print("END!")



