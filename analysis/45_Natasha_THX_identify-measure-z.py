import pandas as pd
import skimage.io as skio
import shared.dataframe as dat
import shared.image as img
import numpy as np
import napari

# input parameters
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220325_Natasha_THZ1/"
sample = 'DMSO'
pixel_size = 102  # nm (sp8 confocal 3144x3144:58.7, mischel 2048x2048:102)
cell_avg_size = 10  # um (Colo)
n_nuclear_convex_dilation = 3
nuclear_centroid_searching_range = 25  # pixel
local_size = 100

# SET UP PARAMETERS
total_fov = 1
total_z = 22

# LOAD IMAGE
im_stack_nuclear = skio.imread("%s%s/%s_RAW_ch00_fov3.tif" % (master_folder, sample, sample), plugin="tifffile")
im_stack_DNAFISH = skio.imread("%s%s/%s_RAW_ch01_fov3.tif" % (master_folder, sample, sample), plugin="tifffile")
im_stack_seg = skio.imread("%s%s/%s_seg_fov3.tif" % (master_folder, sample, sample), plugin="tifffile")

data_z = pd.read_csv('%s%s/%s_fov3_intensity.txt' % (master_folder, sample, sample), na_values=['.'], sep='\t')
data_z['centroid_nuclear'] = [dat.str_to_float(data_z['centroid_nuclear'][i]) for i in range(len(data_z))]
data_z['centroid_x'] = [data_z['centroid_nuclear'][i][0] for i in range(len(data_z))]
data_z['centroid_y'] = [data_z['centroid_nuclear'][i][1] for i in range(len(data_z))]

fov = 0

data = pd.DataFrame(columns=['nuclear',
                             'FOV',
                             'z',
                             'label_nuclear',
                             'centroid_nuclear'])

for fov in range(total_fov):
    print("Analyzing FOV %s/%s" % (fov+1, total_fov))
    data_z_fov = data_z[data_z['FOV'] == fov]
    z_analyze = int(total_z/2)
    data_z_fov_analyze = data_z_fov[data_z_fov['z'] == z_analyze]
    n_nuclear = len(data_z_fov[data_z_fov['z'] == z_analyze])

    # identify each nucleus
    for i in range(n_nuclear):
        print("Analyzing nucleus %s/%s" % (i+1, n_nuclear))
        data_temp = data_z_fov[(data_z_fov['centroid_x'] > data_z_fov_analyze['centroid_x'].tolist()[i] -
                                nuclear_centroid_searching_range)
                               & (data_z_fov['centroid_x'] < data_z_fov_analyze['centroid_x'].tolist()[i] +
                                  nuclear_centroid_searching_range)
                               & (data_z_fov['centroid_y'] > data_z_fov_analyze['centroid_y'].tolist()[i] -
                                  nuclear_centroid_searching_range)
                               & (data_z_fov['centroid_y'] < data_z_fov_analyze['centroid_y'].tolist()[i] +
                                  nuclear_centroid_searching_range)]
        if (len(data_temp) == len(set(data_temp['z'].tolist()))) & (len(data_temp) >= 3):
            limit = int(data_temp['mean_intensity_MYC_DNAFISH_in_nucleus'].tolist()[0] * 2)
            mean_int = []
            for j in range(len(data_temp)):
                z = data_temp['z'].tolist()[j]
                label_nuclear = data_temp['label_nuclear'].tolist()[j]
                original_centroid_nuclear = data_temp['centroid_nuclear'].tolist()[j]

                img_nuclear_seg_convex = im_stack_seg[z]
                img_nuclear = im_stack_nuclear[z]
                img_DNAFISH = im_stack_DNAFISH[z]

                position = img.img_local_position(img_nuclear_seg_convex, original_centroid_nuclear, local_size)
                local_nuclear_seg_convex = img.img_local_seg(img_nuclear_seg_convex, position, label_nuclear)
                local_nuclear = img_nuclear.copy()
                local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
                local_DNAFISH = img_DNAFISH.copy()
                local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]
                local_DNAFISH[local_nuclear_seg_convex == 0] = 0
                local_DNAFISH[local_DNAFISH <= limit] = 0

                temp = local_DNAFISH.ravel()
                temp = [i for i in temp if i != 0]
                # print(np.mean(temp))
                mean_int.append(np.mean(temp))
            data_temp['mean_intensity_over_limit'] = mean_int
                # viewer = napari.view_image(local_DNAFISH)
                # napari.run()
            print(data_temp)
            data_mid = data_temp
            data_mid = data_mid.sort_values(by='mean_intensity_over_limit', ascending=False)
            data_mid.reset_index(drop=True, inplace=True)
            data.loc[len(data.index)] = [i, fov, data_mid['z'].tolist()[0], data_mid['label_nuclear'].tolist()[0],
                                         data_mid['centroid_nuclear'].tolist()[0]]

data.to_csv('%s%s/%s_fov3_z.txt' % (master_folder, sample, sample), index=False, sep='\t')

print("DONE!")