import matplotlib.pyplot as plt
from skimage.measure import label, regionprops
import shared.segmentation as seg
from skimage.morphology import binary_dilation, binary_erosion, disk, dilation, medial_axis
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
import shared.dataframe as dat

# input parameters
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220420_sp8_DMandBRD4_plate/DM_6z_maxProject_imageJ/"
prefix = 'DM_1um-6z-5um_6pos_max_projection_RAW'
sample = 'DM'
pixel_size = 58.7  # nm (sp8 confocal 3144x3144)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
solidity_threshold_nuclear = 0.9
n_nuclear_convex_dilation = 3

# LOAD IMAGE
im_stack_nuclear = skio.imread("%s%s_ch00.tif" % (master_folder, prefix), plugin="tifffile")
im_stack_DNAFISH = skio.imread("%s%s_ch02.tif" % (master_folder, prefix), plugin="tifffile")
im_stack_IF = skio.imread("%s%s_ch01.tif" % (master_folder, prefix), plugin="tifffile")

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
k_dots = 20000  # need to optimize
# segmentation
local_factor_nuclear = 99  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
extreme_val_ecDNA = 18000  # need to be determined
connecting_factor_ecDNA = 10
# radial distribution analysis
radial_interval = 1
radial_max = int(0.8 * local_size)  # ~120
relative_radial_interval = 0.01

# CREATE data
data = pd.DataFrame(columns=['FOV',
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
                             'number_ecDNA',
                             'area_individual_ecDNA',
                             'total_area_ecDNA',
                             'area_percentage_ecDNA_in_nucleus',
                             'mean_area_ecDNA',
                             'max_area_ecDNA',
                             'mean_intensity_individual_ecDNA',
                             'total_intensity_individual_ecDNA',
                             'total_intensity_ecDNA',
                             'participating_coefficient_ecDNA',
                             'centroid_individual_ecDNA',
                             'localization_from_nuclear_centroid_individual_ecDNA',
                             'distance_from_nuclear_centroid_individual_ecDNA',
                             'distance_from_nuclear_edge_individual_ecDNA',
                             'mean_intensity_MYC_IF',
                             'total_intensity_MYC_IF',
                             'g',
                             'dg',
                             'g_value'
                             'radial_distribution_from_centroid',
                             'radial_distribution_from_edge',
                             'radial_distribution_relative_r',
                             'R35l'])

# IMAGING ANALYSIS
for fov in range(total_fov):
    print("Start analyzing FOV %s/%s" % (fov+1, total_fov))

    img_nuclear = im_stack_nuclear[fov]
    img_DNAFISH = im_stack_DNAFISH[fov]
    img_IF = im_stack_IF[fov]

    # nuclear segmentation
    print("Start segmentation...")

    img_nuclear_seg = seg.nuclear_seg1(img_nuclear, local_factor=local_factor_nuclear, min_size=min_size_nuclear,
                                       max_size=max_size_nuclear)
    img_nuclear_seg = seg.filter_solidity(img_nuclear_seg, threshold=solidity_threshold_nuclear)
    img_nuclear_seg_convex = seg.obj_to_convex(img_nuclear_seg)
    img_nuclear_seg_convex = dilation(img_nuclear_seg_convex, disk(n_nuclear_convex_dilation))

    # ecDNA segmentation
    _, img_DNAFISH_seg = seg.find_organelle(img_DNAFISH, 'na', extreme_val=extreme_val_ecDNA,
                                            bg_val=seg.get_bg_int([img_DNAFISH])[0], min_size=10,
                                            max_size=max_size_nuclear)  # min_size: at least needs to have two pixels
    # connect neighbour ecDNA signal
    img_DNAFISH_seg = binary_dilation(img_DNAFISH_seg, disk(connecting_factor_ecDNA))
    for i in range(connecting_factor_ecDNA):
        img_DNAFISH_seg = binary_erosion(img_DNAFISH_seg)
    img_DNAFISH_seg = ndimage.binary_fill_holes(img_DNAFISH_seg)
    img_DNAFISH_seg = obj.remove_small(img_DNAFISH_seg, min_size=20)

    # filter
    img_nuclear_seg_convex = seg.filter_mean_int(img_nuclear_seg_convex, img_DNAFISH_seg, img_DNAFISH, 20000)
    img_DNAFISH_seg[img_nuclear_seg_convex == 0] = 0

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

        ###### HAVEN"T EDIT #####
        # ecDNA related
        FISH_mean_intensity_nuclear = DNAFISH_props[i].mean_intensity
        FISH_total_intensity_nuclear = area_nuclear * FISH_mean_intensity_nuclear
        img_FISH_seg_temp = img_DNAFISH_seg.copy()
        img_FISH_seg_temp[img_nuclear_seg_convex != DNAFISH_props[i].label] = 0
        ecDNA_props = regionprops(label(img_FISH_seg_temp), img_FISH)
        ecDNA_number = len(ecDNA_props)
        ecDNA_area = [ecDNA_props[j].area for j in range(ecDNA_number)]
        ecDNA_total_area = sum(ecDNA_area)
        area_ratio = ecDNA_total_area / area_nuclear
        ecDNA_mean_area = ecDNA_total_area/ecDNA_number
        ecDNA_max_area = max(ecDNA_area)
        ecDNA_mean_int = [ecDNA_props[j].mean_intensity for j in range(ecDNA_number)]
        ecDNA_intensity = [ecDNA_mean_int[j]*ecDNA_area[j] for j in range(ecDNA_number)]
        ecDNA_total_intensity = sum(ecDNA_intensity)
        ecDNA_participating_coefficient = ecDNA_total_intensity * 1.0 / FISH_total_intensity_nuclear
        ecDNA_centroid = [ecDNA_props[j].centroid for j in range(ecDNA_number)]
        ecDNA_localization_from_centroid = [tuple(np.array(ecDNA_centroid[j]) - np.array(centroid_nuclear))
                                            for j in range(ecDNA_number)]
        ecDNA_distance_from_centroid = \
            [(ecDNA_localization_from_centroid[j][0]**2+ecDNA_localization_from_centroid[j][1]**2)**0.5
             for j in range(ecDNA_number)]
        _, distance_map = medial_axis(img_nuclear_seg_convex, return_distance=True)
        ecDNA_distance_from_edge = [distance_map[round(ecDNA_centroid[j][0])][round(ecDNA_centroid[j][1])]
                                    for j in range(ecDNA_number)]

        # MYC expression
        MYC_mean_intensity = IF_props[i].mean_intensity
        MYC_total_intensity = area_nuclear * MYC_mean_intensity

        # auto-correlation
        position = img.img_local_position(img_nuclear_seg_convex, centroid_nuclear, local_size)
        nuclear_seg = img.img_local_seg(img_nuclear_seg_convex, position, i + 1)
        img_FISH_temp = img_FISH.copy()
        FISH = img_FISH_temp[position[0]:position[1], position[2]:position[3]]
        # plt.imsave('%sFISH_fov%s_i%s.tiff' % (master_folder, fov, i), FISH)
        vector = []
        vector_cum_weight = []
        weight = 0
        for m in range(len(nuclear_seg)):
            for n in range(len(nuclear_seg[0])):
                if nuclear_seg[m][n] == 1:
                    vector.append([m, n])
                    if FISH[m][n] > int_thresh_auto_correlation:
                        weight = weight + FISH[m][n] - int_thresh_auto_correlation
                    vector_cum_weight.append(weight)
        random_dot = random.choices(vector, cum_weights=vector_cum_weight, k=100)
        img_dot = np.zeros_like(nuclear_seg)
        for j in random_dot:
            img_dot[j[0]][j[1]] = img_dot[j[0]][j[1]] + 1

        _, r, g, dg = mat.auto_correlation(img_dot, nuclear_seg, rmax)
        # _, r, g, dg = mat.auto_correlation(FISH, nuclear_seg, rmax)

        # radial distribution
        local_nuclear_props = regionprops(label(nuclear_seg))
        local_nuclear_centroid = local_nuclear_props[0].centroid
        pixel = pd.DataFrame(columns=['pixel_distance_from_centroid', 'pixel_distance_from_edge', 'pixel_relative_r',
                                      'pixel_FISH_intensity'])
        _, local_edge_distance_map = medial_axis(nuclear_seg, return_distance=True)
        local_centroid_distance_map = img.distance_map_from_point(nuclear_seg, local_nuclear_centroid)
        local_centroid_distance_map[nuclear_seg == 0] = 0
        local_edge_distance_map[nuclear_seg == 0] = -1
        local_relative_r_map = local_centroid_distance_map/(local_centroid_distance_map+local_edge_distance_map)

        radial_distribution_from_centroid = \
            img.radial_distribution_from_distance_map(nuclear_seg, local_centroid_distance_map, FISH, radial_interval,
                                                      radial_max)
        radial_distribution_from_edge = \
            img.radial_distribution_from_distance_map(nuclear_seg, local_edge_distance_map, FISH, radial_interval,
                                                      radial_max)
        radial_distribution_relative_r = \
            img.radial_distribution_from_distance_map(nuclear_seg, local_relative_r_map, FISH, relative_radial_interval,
                                                      1)
        R35l = img.radial_percentage_from_distance_map(nuclear_seg, local_relative_r_map, FISH, [0, 0.35])
        R35r = img.radial_percentage_from_distance_map(nuclear_seg, local_relative_r_map, FISH, [0.35, 1])
        if R35l > 0.18:
            print("fov%s/i%s" % (fov, i))
            viewer = napari.Viewer()
            viewer.add_image(FISH)
            napari.run()

        # average plot and intensity saturation measurement
        FISH_seg = FISH.copy()
        FISH_seg[nuclear_seg == 0] = 0
        mean_intensity_FISH = np.sum(FISH_seg)/np.sum(nuclear_seg)
        direction = list((np.array(center)-np.array(local_nuclear_centroid)).astype(int))
        FISH_sum = img.sum_up_image(FISH_sum, FISH_seg, direction, 5000.0/mean_intensity_FISH)
        print('max FISH intensity: %s' % np.max(FISH_seg))
        temp = nuclear_seg.copy()
        temp[FISH_seg < 65535] = 0
        print('pixels that are saturated/ FISH: %s' % np.sum(temp))
        print('mean FISH intensity: %s' % np.mean(FISH_seg))

        img_nuclear_temp = img_nuclear.copy()
        local_nuclear = img_nuclear_temp[position[0]:position[1], position[2]:position[3]]
        nuclear_seg_temp = img_nuclear_seg_convex[position[0]:position[1], position[2]:position[3]]
        local_nuclear[nuclear_seg == 0] = 0
        mean_intensity_nuclear = np.sum(local_nuclear)/np.sum(nuclear_seg)
        nuclear_sum = img.sum_up_image(nuclear_sum, local_nuclear, direction, 5000.0/mean_intensity_nuclear)
        print('max nuclear intensity: %s' % np.max(local_nuclear))
        temp = nuclear_seg.copy()
        temp[local_nuclear < 65535] = 0
        print('pixels that are saturated/ nuclear: %s' % np.sum(temp))
        print('mean nuclear intensity: %s' % np.mean(local_nuclear))
        # nuclear_scaled_temp = np.array(nuclear_sum * 65535 / math.ceil(np.max(nuclear_sum))).astype(int)
        # plt.imsave('%savg_nuclear_temp_%s_fov%s_i%s_int.tiff' % (master_folder, sample, fov, i), local_nuclear)
        # plt.imsave('%savg_nuclear_temp_%s_fov%s_i%s.tiff' % (master_folder, sample, fov, i), nuclear_scaled_temp)

        nuclear_seg_sum = img.sum_up_image(nuclear_seg_sum, nuclear_seg, direction, 1)

        img_MYC_temp = img_MYC.copy()
        local_MYC = img_MYC_temp[position[0]:position[1], position[2]:position[3]]
        local_MYC[nuclear_seg == 0] = 0
        print('max MYC intensity: %s' % np.max(local_MYC))
        temp = nuclear_seg.copy()
        temp[local_MYC < 65535] = 0
        print('pixels that are saturated/ MYC: %s' % np.sum(temp))
        print('mean MYC intensity: %s' % np.mean(local_MYC))

        data.loc[len(data.index)] = \
            [fov, label_nuclear, centroid_nuclear, area_nuclear, major_axis_nuclear, minor_axis_nuclear,
             major_and_minor_axis_ratio_nuclear, circularity_nuclear, eccentricity_nuclear, mean_intensity_hoechst,
             total_intensity_hoechst, FISH_mean_intensity_nuclear,
             FISH_total_intensity_nuclear, ecDNA_number, ecDNA_area, ecDNA_total_area, area_ratio, ecDNA_mean_area,
             ecDNA_max_area, ecDNA_mean_int,
             ecDNA_intensity, ecDNA_total_intensity, ecDNA_participating_coefficient, ecDNA_centroid,
             ecDNA_localization_from_centroid, ecDNA_distance_from_centroid, ecDNA_distance_from_edge,
             MYC_mean_intensity, MYC_total_intensity,
             g, dg, radial_distribution_from_centroid, radial_distribution_from_edge, radial_distribution_relative_r,
             R35l, R35r]

viewer = napari.Viewer()
viewer.add_image(img_nuclear)
napari.run()
viewer = napari.Viewer()
viewer.add_image(img_FISH)
napari.run()
viewer = napari.Viewer()
viewer.add_image(img_MYC)
napari.run()

nuclear_average_mean_intensity = data['nuclear_mean_intensity'].mean()
FISH_average_mean_intensity_nuclear = data['FISH_mean_intensity_nuclear'].mean()
FISH_scaled = np.array(FISH_sum * FISH_average_mean_intensity_nuclear / 5000.0).astype(int)
nuclear_scaled = np.array(nuclear_sum * nuclear_average_mean_intensity / 5000.0).astype(int)
nuclear_seg_scaled = np.array(nuclear_seg_sum*65535/math.ceil(np.max(nuclear_seg_sum))).astype(int)
plt.imsave('%savg_FISH_%s.tiff' % (master_folder, sample), FISH_scaled)
plt.imsave('%savg_nuclear_%s.tiff' % (master_folder, sample), nuclear_scaled)
plt.imsave('%savg_nuclear_seg_%s.tiff' % (master_folder, sample), nuclear_seg_scaled)

data.to_csv('%s%s.txt' % (master_folder, sample), index=False, sep='\t')

print("DONE!")




