import matplotlib.pyplot as plt
from skimage.measure import label, regionprops
from skimage.morphology import medial_axis, dilation, erosion
import pandas as pd
import math
import numpy as np
import shared.image as img
import random
import shared.math as mat
import napari
import skimage.io as skio
import shared.dataframe as dat
import shared.objects as obj

# input parameters
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220301_ecDNA_ctrlAndJQ1_NatashaFile/"
prefix = '210921_COLODM_washout_mycFISH'
total_fov = 2
sample = 'DMSO'
pixel_size = 40  # nm (Zeiss confocal scope)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
nuclear_centroid_searching_range = 50  # pixel
z_analyze = [16, 8]
# 3hrJQ1_3hr1uMtriptolide_WO [22, 15]
# 3hrJQ1_3hrDMSO_WO [25, 25, 16]
# DMSO3hr [9, 14, 9]
# JQ13hr [8, 9]

# DMSO [16, 8]
# JQ13hr [21, 19, 18]

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
int_thresh_auto_correlation = 300  # need to be determined
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
                             'area_nuclear',
                             'major_axis_nuclear',
                             'minor_axis_nuclear',
                             'major_and_minor_axis_ratio_nuclear',
                             'circularity_nuclear',
                             'eccentricity_nuclear',
                             'mean_intensity_hoechst',
                             'total_intensity_hoechst',
                             'mean_intensity_MYC_DNAFISH_in_nucleus',
                             'total_intensity_MYC_DNAFISH_in_nucleus',
                             'g',
                             'g_correct',
                             'dg',
                             'g_value',
                             'g_correct_value',
                             'radial_distribution_from_centroid',
                             'radial_distribution_from_edge',
                             'radial_distribution_relative_r',
                             'R35l'])

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
                for j in range(1):
                    print("Analyzing sub-nucleus %s/3" % (j+1))
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
                    plt.imsave('%sFISH_fov%s_i%s_j%s_z%s.tiff' % (master_folder, fov, i, j, z_nuclear_temp), local_DNAFISH)

                    # auto correlation
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
                    """img_dot_remove_bg = dilation(img_dot)
                    img_dot_remove_bg = dilation(img_dot_remove_bg)
                    img_dot_remove_bg = erosion(img_dot_remove_bg)
                    img_dot_remove_bg = erosion(img_dot_remove_bg)
                    img_dot_remove_bg = erosion(img_dot_remove_bg)
                    img_dot_remove_bg = dilation(img_dot_remove_bg)
                    img_dot_remove_bg_seg_filter = obj.remove_small(label(img_dot_remove_bg), 35)
                    img_dot_remove_bg_filter = img_dot_remove_bg.copy()
                    img_dot_remove_bg_filter[img_dot_remove_bg_seg_filter == 0] = 0"""
                    """viewer = napari.Viewer()
                    viewer.add_image(img_dot)
                    napari.run()"""

                    _, r, g, dg = mat.auto_correlation(img_dot, local_nuclear_seg_convex, rmax)
                    g_value = (g[1] + g[2] + g[3] + g[4] + g[5]) * 0.2

                    local_DNAFISH_temp = local_DNAFISH.copy()
                    local_DNAFISH_temp[local_nuclear_seg_convex == 0] = 0
                    total_intensity_MYC_DNAFISH_in_dot = np.sum(local_DNAFISH_temp)
                    g_correct = list(np.array(g) * (total_intensity_MYC_DNAFISH_in_dot / 10000000.0))
                    g_correct_value = g_value * (total_intensity_MYC_DNAFISH_in_dot / 10000000.0)

                    plt.imsave('%simg_dot_fov%s_i%s_j%s_z%s_g%s_g_correct%s.tiff' %
                               (master_folder, fov, i, j, z_nuclear_temp, g_value, g_correct_value), img_dot)
                    """plt.imsave('%simg_dot_remove_bg_filter_fov%s_i%s_j%s_z%s_g%s_g_correct%s.tiff' %
                               (master_folder, fov, i, j, z_nuclear_temp, g_value, g_correct_value), img_dot_remove_bg_filter)"""

                    """viewer = napari.Viewer()
                    viewer.add_image(local_DNAFISH)
                    napari.run()
                    viewer = napari.Viewer()
                    viewer.add_image(img_dot)
                    napari.run()"""

                    # radial distribution
                    local_nuclear_props = regionprops(label(local_nuclear_seg_convex))
                    local_nuclear_centroid = local_nuclear_props[0].centroid
                    _, local_edge_distance_map = medial_axis(local_nuclear_seg_convex, return_distance=True)
                    local_centroid_distance_map = img.distance_map_from_point(local_nuclear_seg_convex,
                                                                              local_nuclear_centroid)
                    local_centroid_distance_map[local_nuclear_seg_convex == 0] = 0
                    local_edge_distance_map[local_nuclear_seg_convex == 0] = -1
                    local_relative_r_map = local_centroid_distance_map / (
                                local_centroid_distance_map + local_edge_distance_map)

                    radial_distribution_from_centroid = \
                        img.radial_distribution_from_distance_map(local_nuclear_seg_convex, local_centroid_distance_map,
                                                                  local_DNAFISH, radial_interval, radial_max)
                    radial_distribution_from_edge = \
                        img.radial_distribution_from_distance_map(local_nuclear_seg_convex, local_edge_distance_map,
                                                                  local_DNAFISH, radial_interval, radial_max)
                    radial_distribution_relative_r = \
                        img.radial_distribution_from_distance_map(local_nuclear_seg_convex, local_relative_r_map,
                                                                  local_DNAFISH, relative_radial_interval, 1)
                    R35l = img.radial_percentage_from_distance_map(local_nuclear_seg_convex, local_relative_r_map,
                                                                   local_DNAFISH, [0, 0.35])

                    # average plot and intensity saturation measurement
                    local_DNAFISH_seg = local_DNAFISH.copy()
                    local_DNAFISH_seg[local_nuclear_seg_convex == 0] = 0
                    mean_intensity_DNAFISH = np.sum(local_DNAFISH_seg) / np.sum(local_nuclear_seg_convex)
                    direction = list((np.array(avg_img_center) - np.array(local_nuclear_centroid)).astype(int))
                    avg_img_DNAFISH = img.sum_up_image(avg_img_DNAFISH, local_DNAFISH_seg, direction,
                                                       500.0 / mean_intensity_DNAFISH)
                    print('max FISH intensity: %s' % np.max(local_DNAFISH_seg))
                    temp = local_nuclear_seg_convex.copy()
                    temp[local_DNAFISH_seg < 65535] = 0
                    print('pixels that are saturated/ FISH: %s' % np.sum(temp))
                    print('mean FISH intensity: %s' % np.mean(local_DNAFISH_seg))

                    local_nuclear_seg = local_nuclear.copy()
                    local_nuclear_seg[local_nuclear_seg_convex == 0] = 0
                    mean_intensity_nuclear = np.mean(local_nuclear_seg)
                    avg_img_nuclear = img.sum_up_image(avg_img_nuclear, local_nuclear_seg, direction,
                                                       500.0 / mean_intensity_nuclear)
                    print('max nuclear intensity: %s' % np.max(local_nuclear))
                    temp = local_nuclear_seg_convex.copy()
                    temp[local_nuclear < 65535] = 0
                    print('pixels that are saturated/ nuclear: %s' % np.sum(temp))
                    print('mean nuclear intensity: %s' % np.mean(local_nuclear))

                    avg_img_nuclear_seg = img.sum_up_image(avg_img_nuclear_seg, local_nuclear_seg_convex, direction, 1)

                    # nuclear morphology
                    area_nuclear = data_temp['area_nuclear'][j]
                    major_axis_nuclear = local_nuclear_props[0].major_axis_length
                    minor_axis_nuclear = local_nuclear_props[0].minor_axis_length
                    major_and_minor_axis_ratio_nuclear = major_axis_nuclear*1.0/minor_axis_nuclear
                    circularity_nuclear = (4 * math.pi * area_nuclear) / (local_nuclear_props[0].perimeter ** 2)
                    eccentricity_nuclear = local_nuclear_props[0].eccentricity

                    data.loc[len(data.index)] = [fov, z_analyze[fov], i, z_nuclear_temp, label_nuclear_temp,
                                                 original_centroid_nuclear_temp,
                                                 area_nuclear, major_axis_nuclear, minor_axis_nuclear,
                                                 major_and_minor_axis_ratio_nuclear, circularity_nuclear,
                                                 eccentricity_nuclear, mean_intensity_nuclear,
                                                 mean_intensity_nuclear * area_nuclear, mean_intensity_DNAFISH,
                                                 mean_intensity_DNAFISH * area_nuclear, g, g_correct, dg, g_value,
                                                 g_correct_value,
                                                 radial_distribution_from_centroid, radial_distribution_from_edge,
                                                 radial_distribution_relative_r, R35l]

data.to_csv('%s%s.txt' % (master_folder, sample), index=False, sep='\t')

nuclear_average_mean_intensity = data['mean_intensity_hoechst'].mean()
DNAFISH_average_mean_intensity_nuclear = data['mean_intensity_MYC_DNAFISH_in_nucleus'].mean()
DNAFISH_scaled = np.array(avg_img_DNAFISH * DNAFISH_average_mean_intensity_nuclear / 500.0).astype(int)
nuclear_scaled = np.array(avg_img_nuclear * nuclear_average_mean_intensity / 500.0).astype(int)
nuclear_seg_scaled = np.array(avg_img_nuclear_seg * 65535/math.ceil(np.max(avg_img_nuclear_seg))).astype(int)
plt.imsave('%savg_DNAFISH_%s.tiff' % (master_folder, sample), DNAFISH_scaled)
plt.imsave('%savg_nuclear_%s.tiff' % (master_folder, sample), nuclear_scaled)
plt.imsave('%savg_nuclear_seg_%s.tiff' % (master_folder, sample), nuclear_seg_scaled)

"""data_mean = pd.DataFrame(columns=['FOV',
                                  'z_analyze',
                                  'nuclear',
                                  'centroid_nuclear',
                                  'area_nuclear',
                                  'major_axis_nuclear',
                                  'minor_axis_nuclear',
                                  'major_and_minor_axis_ratio_nuclear',
                                  'circularity_nuclear',
                                  'eccentricity_nuclear',
                                  'mean_intensity_hoechst',
                                  'total_intensity_hoechst',
                                  'mean_intensity_MYC_DNAFISH_in_nucleus',
                                  'total_intensity_MYC_DNAFISH_in_nucleus',
                                  'g',
                                  'g_value',
                                  'radial_distribution_relative_r',
                                  'R35l'])

for x in range(int(len(data)/3)):
    g_mean = list((np.array(data['g'].tolist()[3*x]) + np.array(data['g'].tolist()[3*x+1])
                   + np.array(data['g'].tolist()[3*x+2]))*1.0/3)
    g_value_mean = (data['g_value'][3*x] + data['g_value'][3*x+1] + data['g_value'][3*x+2]) *1.0/3
    radial_distribution_relative_r_mean = list((np.array(data['radial_distribution_relative_r'].tolist()[3*x])
                                                + np.array(data['radial_distribution_relative_r'].tolist()[3*x+1])
                                                + np.array(data['radial_distribution_relative_r'].tolist()[3*x+2]))*1.0/3)
    R35l_mean = (data['R35l'][3*x] + data['R35l'][3*x+1] + data['R35l'][3*x+2]) *1.0/3
    data_mean.loc[len(data_mean.index)] = [data['FOV'][3*x], data['z_analyze'][3*x],
                                           data['nuclear'][3*x], data['centroid_nuclear'][3*x],
                                           data['area_nuclear'][3*x], data['major_axis_nuclear'][3*x],
                                           data['minor_axis_nuclear'][3*x],
                                           data['major_and_minor_axis_ratio_nuclear'][3*x],
                                           data['circularity_nuclear'][3*x], data['eccentricity_nuclear'][3*x],
                                           data['mean_intensity_hoechst'][3*x], data['total_intensity_hoechst'][3*x],
                                           data['mean_intensity_MYC_DNAFISH_in_nucleus'][3*x],
                                           data['total_intensity_MYC_DNAFISH_in_nucleus'][3*x], g_mean, g_value_mean,
                                           radial_distribution_relative_r_mean, R35l_mean]

data_mean.to_csv('%s%s_mean.txt' % (master_folder, sample), index=False, sep='\t')"""

"""viewer = napari.Viewer()
viewer.add_image(im_stack)
napari.run()
viewer = napari.Viewer()
viewer.add_image(nuclear_seg_convex)
napari.run()"""

print("DONE!")
