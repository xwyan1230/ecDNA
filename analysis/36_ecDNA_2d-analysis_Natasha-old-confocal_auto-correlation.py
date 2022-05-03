import matplotlib.pyplot as plt
from skimage.measure import label, regionprops
from skimage.morphology import medial_axis
import pandas as pd
import math
import numpy as np
import shared.image as img
import random
import shared.math as mat
import napari
import skimage.io as skio
import shared.dataframe as dat

# input parameters
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/211102_3hr_JQ1washout/JQ13hr/"
prefix = '211102_COLODM'
total_fov = 2
sample = 'JQ13hr'
pixel_size = 40  # nm (Zeiss confocal scope)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
nuclear_centroid_searching_range = 50  # pixel
z_analyze = [8, 9]
# 3hrJQ1_3hr1uMtriptolide_WO [22, 15]
# 3hrJQ1_3hrDMSO_WO [25, 25, 16]
# DMSO3hr [9, 14, 9]
# JQ13hr [8, 9]

# SET UP PARAMETERS
# single nuclear analysis
local_size = int(0.9 * cell_avg_size * 1000/pixel_size)  # single FOV only, assuming squared pixel, ~150
# avg_img
avg_img_size = 2 * local_size  # for generating avg_img, ~300
avg_img_DNAFISH = np.zeros(shape=(avg_img_size, avg_img_size))
avg_img_nuclear = np.zeros(shape=(avg_img_size, avg_img_size))
avg_img_nuclear_seg = np.zeros(shape=(avg_img_size, avg_img_size))
avg_img_center = [int(avg_img_size / 2), int(avg_img_size / 2)]
# auto-correlation analysis
rmax = int(0.67 * local_size)  # ~100
int_thresh_auto_correlation = 500  # need to be determined
k_dots = 10000  # need to optimize
# segmentation
local_factor_nuclear = rmax if (rmax % 2 == 1) else rmax+1  # ~99, needs to be odd number
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
extreme_val_ecDNA = 18000  # need to be determined
connecting_factor_ecDNA = 10
# radial distribution analysis
radial_interval = 1
radial_max = int(0.8 * local_size)  # ~120
relative_radial_interval = 0.01

data = pd.DataFrame(columns=['FOV',
                             'z_analyze',
                             'nuclear',
                             'z',
                             'label_nuclear',
                             'centroid_nuclear',
                             'int_thresh_auto_correlation',
                             'g',
                             'dg',
                             'g_value'])

data_g = pd.DataFrame(columns=['FOV', 'nuclear', 'g_value'])

# OPEN Z FILE
data_z = pd.read_csv('%s%s_intensity.txt' % (master_folder, sample), na_values=['.'], sep='\t')
data_z['centroid_nuclear'] = [dat.str_to_float(data_z['centroid_nuclear'][i]) for i in range(len(data_z))]
data_z['centroid_x'] = [data_z['centroid_nuclear'][i][0] for i in range(len(data_z))]
data_z['centroid_y'] = [data_z['centroid_nuclear'][i][1] for i in range(len(data_z))]

# IMAGING ANALYSIS
for fov in range(total_fov):
    print("Analyzing FOV %s/%s" % (fov+1, total_fov))
    data_z_fov = data_z[data_z['FOV'] == fov]
    z_max = 0
    z_cell = 0
    for z in range(max(data_z_fov['z'])):
        temp = len(data_z_fov[data_z_fov['z'] == z])
        if temp > z_cell:
            z_cell = temp
            z_max = z
    print(z_max)
    data_z_fov_max = data_z_fov[data_z_fov['z'] == z_analyze[fov]]
    # load images
    im_stack = skio.imread("%s%s_%s_%s.tif" % (master_folder, prefix, sample, fov + 1), plugin="tifffile")
    nuclear_seg_convex = skio.imread('%snuclear_seg_convex_%s_%s.tif' % (master_folder, sample, fov + 1),
                                     plugin="tifffile")
    # identify each nucleus
    for i in range(len(data_z_fov_max)):
        print("Analyzing nucleus %s/%s" % (i+1, len(data_z_fov_max)))
        data_temp = data_z_fov[(data_z_fov['centroid_x'] > data_z_fov_max['centroid_x'].tolist()[i] - nuclear_centroid_searching_range) &
                               (data_z_fov['centroid_x'] < data_z_fov_max['centroid_x'].tolist()[i] + nuclear_centroid_searching_range) &
                               (data_z_fov['centroid_y'] > data_z_fov_max['centroid_y'].tolist()[i] - nuclear_centroid_searching_range) &
                               (data_z_fov['centroid_y'] < data_z_fov_max['centroid_y'].tolist()[i] + nuclear_centroid_searching_range)]
        if (len(data_temp) == len(set(data_temp['z'].tolist()))) & (len(data_temp) >= 3):
            data_temp = data_temp.sort_values(by='total_intensity_MYC_DNAFISH_in_nucleus', ascending=False)
            data_temp.reset_index(drop=True, inplace=True)
            if max(data_temp[0:3]['z']) - min(data_temp[0:3]['z']) == 2:
                j = 0
                z_nuclear_temp = data_temp['z'][j]
                label_nuclear_temp = data_temp['label_nuclear'][j]
                original_centroid_nuclear_temp = data_temp['centroid_nuclear'][j]

                img_nuclear_seg_convex = nuclear_seg_convex[z_nuclear_temp]
                img_nuclear = im_stack[z_nuclear_temp][1]
                img_DNAFISH = im_stack[z_nuclear_temp][0]

                position = img.img_local_position(img_nuclear_seg_convex, original_centroid_nuclear_temp, local_size)
                local_nuclear_seg_convex = img.img_local_seg(img_nuclear_seg_convex, position, label_nuclear_temp)
                local_nuclear = img_nuclear.copy()
                local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
                local_DNAFISH = img_DNAFISH.copy()
                local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]
                plt.imsave('%sFISH_fov%s_i%s_z%s.tiff' % (master_folder, fov, i, z_nuclear_temp), local_DNAFISH)
                plt.imsave('%slocal_seg_fov%s_i%s_z%s.tiff' % (master_folder, fov, i, z_nuclear_temp),
                           local_nuclear_seg_convex)

                # auto correlation
                g_list = []
                for int_thresh_auto_correlation in np.arange(0, 2000, 100):
                    vector = []
                    vector_cum_weight = []
                    weight = 0
                    for m in range(len(local_nuclear_seg_convex)):
                        for n in range(len(local_nuclear_seg_convex[0])):
                            if local_nuclear_seg_convex[m][n] == 1:
                                vector.append([m, n])
                                if local_DNAFISH[m][n] > int_thresh_auto_correlation:
                                    weight = weight + local_DNAFISH[m][n] - int_thresh_auto_correlation
                                vector_cum_weight.append(weight)
                    random_dot = random.choices(vector, cum_weights=vector_cum_weight, k=k_dots)
                    img_dot = np.zeros_like(local_nuclear_seg_convex)
                    for l in random_dot:
                        img_dot[l[0]][l[1]] = img_dot[l[0]][l[1]] + 1

                    _, r, g, dg = mat.auto_correlation(img_dot, local_nuclear_seg_convex, rmax)
                    g_value = (g[1] + g[2] + g[3] + g[4] + g[5] + g[6] + g[7] + g[8]) * 1.0 / 8
                    g_list.append(g_value)
                    plt.imsave('%simg_dot_fov%s_i%s_z%s_thresh%s_g%s.tiff' %
                               (master_folder, fov, i, z_nuclear_temp, int_thresh_auto_correlation, g_value), img_dot)

                    data.loc[len(data.index)] = [fov, z_analyze[fov], i, z_nuclear_temp, label_nuclear_temp,
                                                 original_centroid_nuclear_temp, int_thresh_auto_correlation, g, dg,
                                                 g_value]
                data_g.loc[len(data_g.index)] = [fov, i, g_list]

data.to_csv('%s%s_different_int_threshold.txt' % (master_folder, sample), index=False, sep='\t')
data_g.to_csv('%s%s_different_int_threshold_g.txt' % (master_folder, sample), index=False, sep='\t')

print("DONE!")
