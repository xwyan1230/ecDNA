import matplotlib.pyplot as plt
from skimage.measure import label, regionprops
import shared.segmentation as seg
from skimage.morphology import binary_dilation, binary_erosion, disk, dilation, medial_axis, erosion
import shared.objects as obj
from scipy import ndimage
import pandas as pd
import math
import numpy as np
import shared.image as img
import random
import shared.math as mat
import napari
import skimage.io as skio
import tifffile as tif
import shared.dataframe as dat

# input parameters

"""master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220407_sp8_DMandHSR/HSR_singleZ/"
# prefix = 'DM_singleZ_25pos_RAW'
sample = 'HSR'"""
"""master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220420_sp8_DMandBRD4_plate/DM_singleZ_imageJ/"
prefix = 'DM_singleZ_25pos_RAW'
sample = 'DM'"""
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220420_sp8_DMandBRD4_plate/BRD4KO_singleZ_imageJ/"
prefix = 'DM-BRD4KO_singleZ_8pos_RAW'
sample = 'BRD4KO'
pixel_size = 58.7  # nm (sp8 confocal 3144x3144)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
solidity_threshold_nuclear = 0.9
n_nuclear_convex_dilation = 3
DNAFISH_bg_correct_factor = 8000

# LOAD IMAGE
im_stack_nuclear = skio.imread("%s%s_ch00.tif" % (master_folder, prefix), plugin="tifffile")
im_stack_DNAFISH = skio.imread("%sDNAFISH_corrected_%s_%s.tif" % (master_folder, sample, DNAFISH_bg_correct_factor), plugin="tifffile")
im_stack_IF = skio.imread("%s%s_ch01.tif" % (master_folder, prefix), plugin="tifffile")
im_stack_nuclear_seg_convex = skio.imread("%snuclear_seg_convex_%s.tif" % (master_folder, sample), plugin="tifffile")

# SET UP PARAMETERS
total_fov = im_stack_nuclear.shape[0]
# single nuclear analysis
local_size = 150  # single FOV only, assuming squared pixel, ~150, int(0.9 * cell_avg_size * 1000/pixel_size)
# avg_img
avg_img_size = 300  # for generating avg_img, ~300, 2 * local_size
avg_img_DNAFISH = np.zeros(shape=(avg_img_size, avg_img_size))
avg_img_nuclear = np.zeros(shape=(avg_img_size, avg_img_size))
avg_img_nuclear_seg = np.zeros(shape=(avg_img_size, avg_img_size))
avg_img_center = [int(avg_img_size / 2), int(avg_img_size / 2)]
# auto-correlation analysis
rmax = 100  # int(0.67 * local_size)
int_thresh_auto_correlation = 10000  # need to be determined
k_dots = 10000  # need to optimize
# segmentation
local_factor_nuclear = 99  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
extreme_val_ecDNA = 15000  # need to be determined
connecting_factor_ecDNA = 8
ecDNA_seg_min_size1 = 10
ecDNA_seg_min_size2 = 20
ecDNA_filter_intensity = 20000
# radial distribution analysis
radial_interval = 1
radial_max = int(0.8 * local_size)  # ~120
relative_radial_interval = 0.01

# CREATE data
data = pd.DataFrame(columns=['FOV',
                             'DNAFISH_bg_correct_factor',
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
                             'mean_intensity_MYC_IF',
                             'total_intensity_MYC_IF',
                             'g',
                             'g_correct',
                             'dg',
                             'g_value',
                             'g_correct_value',
                             'radial_distribution_relative_r',
                             'R35l'])

# IMAGING ANALYSIS
nuclear_seg_convex_filtered = np.zeros(shape=(im_stack_nuclear.shape[0], im_stack_nuclear.shape[1],
                                              im_stack_nuclear.shape[2]), dtype=np.uint16)
ecDNA_seg = np.zeros(shape=(im_stack_nuclear.shape[0], im_stack_nuclear.shape[1], im_stack_nuclear.shape[2]),
                     dtype=np.uint16)
for fov in range(total_fov):
    print("Start analyzing FOV %s/%s" % (fov + 1, total_fov))
    img_nuclear = im_stack_nuclear[fov]
    img_DNAFISH = im_stack_DNAFISH[fov]
    img_IF = im_stack_IF[fov]
    img_nuclear_seg_convex = im_stack_nuclear_seg_convex[fov]

    """# ecDNA segmentation
    _, img_DNAFISH_seg = seg.find_organelle(img_DNAFISH, 'na', extreme_val=extreme_val_ecDNA,
                                            bg_val=seg.get_bg_int([img_DNAFISH])[0], min_size=ecDNA_seg_min_size1,
                                            max_size=max_size_nuclear)
    img_DNAFISH_seg = binary_dilation(img_DNAFISH_seg, disk(connecting_factor_ecDNA))
    for i in range(connecting_factor_ecDNA):
        img_DNAFISH_seg = binary_erosion(img_DNAFISH_seg)
    img_DNAFISH_seg = ndimage.binary_fill_holes(img_DNAFISH_seg)
    img_DNAFISH_seg = obj.remove_small(img_DNAFISH_seg, min_size=ecDNA_seg_min_size2)
    ecDNA_seg[fov] = img_DNAFISH_seg

    # filter
    img_nuclear_seg_convex = seg.filter_mean_int(img_nuclear_seg_convex, img_DNAFISH_seg, img_DNAFISH,
                                                 ecDNA_filter_intensity)
    img_DNAFISH_seg[img_nuclear_seg_convex == 0] = 0
    nuclear_seg_convex_filtered[fov] = img_nuclear_seg_convex"""

    # props
    print("Start analyzing features...")
    DNAFISH_props = regionprops(img_nuclear_seg_convex, img_DNAFISH)
    nuclear_props = regionprops(img_nuclear_seg_convex, img_nuclear)
    IF_props = regionprops(img_nuclear_seg_convex, img_IF)

    for i in range(len(DNAFISH_props)):
        # nuclear morphology
        label_nuclear = DNAFISH_props[i].label
        centroid_nuclear = DNAFISH_props[i].centroid
        area_nuclear = DNAFISH_props[i].area
        major_axis_nuclear = DNAFISH_props[i].major_axis_length
        minor_axis_nuclear = DNAFISH_props[i].minor_axis_length
        major_and_minor_axis_ratio_nuclear = major_axis_nuclear * 1.0 / minor_axis_nuclear
        circularity_nuclear = (4 * math.pi * area_nuclear) / (DNAFISH_props[i].perimeter ** 2)
        eccentricity_nuclear = DNAFISH_props[i].eccentricity
        mean_intensity_hoechst = nuclear_props[i].mean_intensity
        total_intensity_hoechst = area_nuclear * mean_intensity_hoechst

        # ecDNA related
        mean_intensity_MYC_DNAFISH_in_nucleus = DNAFISH_props[i].mean_intensity
        total_intensity_MYC_DNAFISH_in_nucleus = area_nuclear * mean_intensity_MYC_DNAFISH_in_nucleus

        # MYC expression
        mean_intensity_MYC_IF = IF_props[i].mean_intensity
        total_intensity_MYC_IF = area_nuclear * mean_intensity_MYC_IF

        # auto-correlation
        print("Start auto correlation analysis...")
        position = img.img_local_position(img_nuclear_seg_convex, centroid_nuclear, local_size)
        local_nuclear_seg = img.img_local_seg(img_nuclear_seg_convex, position, i + 1)
        img_DNAFISH_temp = img_DNAFISH.copy()
        local_DNAFISH = img_DNAFISH_temp[position[0]:position[1], position[2]:position[3]]
        plt.imsave('%sFISH_fov%s_i%s.tiff' % (master_folder, fov, i), local_DNAFISH)
        vector = []
        vector_cum_weight = []
        weight = 0
        for m in range(len(local_nuclear_seg)):
            for n in range(len(local_nuclear_seg[0])):
                if local_nuclear_seg[m][n] == 1:
                    vector.append([m, n])
                    if local_DNAFISH[m][n] > int_thresh_auto_correlation:
                        weight = weight + local_DNAFISH[m][n] - int_thresh_auto_correlation
                    vector_cum_weight.append(weight)
        random_dot = random.choices(vector, cum_weights=vector_cum_weight, k=k_dots)
        img_dot = np.zeros_like(local_nuclear_seg)
        for j in random_dot:
            img_dot[j[0]][j[1]] = img_dot[j[0]][j[1]] + 1
        """img_dot_remove_bg = dilation(img_dot)
        img_dot_remove_bg = erosion(img_dot_remove_bg)
        img_dot_remove_bg = erosion(img_dot_remove_bg)
        img_dot_remove_bg = dilation(img_dot_remove_bg)"""

        _, r, g, dg = mat.auto_correlation(img_dot, local_nuclear_seg, rmax)
        g_value = (g[1] + g[2] + g[3] + g[4] + g[5]) * 0.2

        local_DNAFISH_temp = local_DNAFISH.copy()
        local_DNAFISH_temp[local_nuclear_seg == 0] = 0
        total_intensity_MYC_DNAFISH_in_dot = np.sum(local_DNAFISH_temp)
        g_correct = list(np.array(g) * (total_intensity_MYC_DNAFISH_in_dot / 10000000.0))
        g_correct_value = g_value * (total_intensity_MYC_DNAFISH_in_dot / 10000000.0)

        plt.imsave('%simg_dot_fov%s_i%s_g%s_g_correct%s.tiff' % (master_folder, fov, i, g_value, g_correct_value), img_dot)
        """plt.imsave('%simg_dot_remove_bg_fov%s_i%s_g%s_g_correct%s.tiff' % (master_folder, fov, i, g_value,
                                                                           g_correct_value), img_dot_remove_bg)"""

        # radial distribution
        print("Start radial distribution analysis...")
        local_nuclear_props = regionprops(label(local_nuclear_seg))
        local_nuclear_centroid = local_nuclear_props[0].centroid
        pixel = pd.DataFrame(columns=['pixel_distance_from_centroid', 'pixel_distance_from_edge', 'pixel_relative_r',
                                      'pixel_FISH_intensity'])
        _, local_edge_distance_map = medial_axis(local_nuclear_seg, return_distance=True)
        local_centroid_distance_map = img.distance_map_from_point(local_nuclear_seg, local_nuclear_centroid)
        local_centroid_distance_map[local_nuclear_seg == 0] = 0
        local_edge_distance_map[local_nuclear_seg == 0] = -1
        local_relative_r_map = local_centroid_distance_map/(local_centroid_distance_map+local_edge_distance_map)

        radial_distribution_relative_r = \
            img.radial_distribution_from_distance_map(local_nuclear_seg, local_relative_r_map, local_DNAFISH,
                                                      relative_radial_interval, 1)
        R35l = img.radial_percentage_from_distance_map(local_nuclear_seg, local_relative_r_map, local_DNAFISH,
                                                       [0, 0.35])

        # average plot and intensity saturation measurement
        print("Start generating average plot...")
        local_DNAFISH_seg = local_DNAFISH.copy()
        local_DNAFISH_seg[local_nuclear_seg == 0] = 0
        mean_intensity_DNAFISH = np.mean(local_DNAFISH_seg)
        direction = list((np.array(avg_img_center)-np.array(local_nuclear_centroid)).astype(int))
        avg_img_DNAFISH = img.sum_up_image(avg_img_DNAFISH, local_DNAFISH_seg, direction, 5000.0/mean_intensity_DNAFISH)
        print('max FISH intensity: %s' % np.max(local_DNAFISH_seg))
        print('mean FISH intensity: %s' % np.mean(local_DNAFISH_seg))

        img_nuclear_temp = img_nuclear.copy()
        local_nuclear = img_nuclear_temp[position[0]:position[1], position[2]:position[3]]
        local_nuclear[local_nuclear_seg == 0] = 0
        mean_intensity_nuclear = np.mean(local_nuclear)
        avg_img_nuclear = img.sum_up_image(avg_img_nuclear, local_nuclear, direction, 5000.0/mean_intensity_nuclear)
        print('max nuclear intensity: %s' % np.max(local_nuclear))
        print('mean nuclear intensity: %s' % np.mean(local_nuclear))

        avg_img_nuclear_seg = img.sum_up_image(avg_img_nuclear_seg, local_nuclear_seg, direction, 1)

        img_IF_temp = img_IF.copy()
        local_IF = img_IF_temp[position[0]:position[1], position[2]:position[3]]
        local_IF[local_nuclear_seg == 0] = 0
        print('max MYC intensity: %s' % np.max(local_IF))
        print('mean MYC intensity: %s' % np.mean(local_IF))

        data.loc[len(data.index)] = \
            [fov, DNAFISH_bg_correct_factor, label_nuclear, centroid_nuclear, area_nuclear, major_axis_nuclear,
             minor_axis_nuclear, major_and_minor_axis_ratio_nuclear, circularity_nuclear, eccentricity_nuclear,
             mean_intensity_hoechst, total_intensity_hoechst, mean_intensity_MYC_DNAFISH_in_nucleus,
             total_intensity_MYC_DNAFISH_in_nucleus, mean_intensity_MYC_IF, total_intensity_MYC_IF, g, g_correct, dg,
             g_value, g_correct_value, radial_distribution_relative_r, R35l]

# tif.imsave('%snuclear_seg_convex_filtered_%s.tif' % (master_folder, sample), nuclear_seg_convex_filtered)
# tif.imsave('%secDNA_seg_%s.tif' % (master_folder, sample), ecDNA_seg)

data.to_csv('%s%s.txt' % (master_folder, sample), index=False, sep='\t')

average_mean_intensity_hoechst = data['mean_intensity_hoechst'].mean()
average_mean_intensity_MYC_DNAFISH_in_nuclear = data['mean_intensity_MYC_DNAFISH_in_nucleus'].mean()
avg_img_DNAFISH_scaled = np.array(avg_img_DNAFISH * average_mean_intensity_MYC_DNAFISH_in_nuclear / 5000.0).astype(int)
avg_img_nuclear_scaled = np.array(avg_img_nuclear * average_mean_intensity_hoechst / 5000.0).astype(int)
avg_img_nuclear_seg_scaled = np.array(avg_img_nuclear_seg*65535/math.ceil(np.max(avg_img_nuclear_seg))).astype(int)
plt.imsave('%savg_DNAFISH_%s.tiff' % (master_folder, sample), avg_img_DNAFISH_scaled)
plt.imsave('%savg_nuclear_%s.tiff' % (master_folder, sample), avg_img_nuclear_scaled)
plt.imsave('%savg_nuclear_seg_%s.tiff' % (master_folder, sample), avg_img_nuclear_seg_scaled)

print("DONE!")




