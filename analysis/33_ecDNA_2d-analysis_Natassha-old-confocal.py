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
import imageio

# input parameters
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/211102_3hr_JQ1washout/"
prefix = '211102_COLODM'
# _3hrJQ1_3hr1uMtriptolide_WO_1
total_fov = 2
sample = '3hrJQ1_3hr1uMtriptolide_WO'
pixel_size = 40  # nm (Zeiss confocal scope)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
solidity_threshold_nuclear = 0.9
n_nuclear_convex_dilation = 3
total_fov = 2

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
int_thresh_auto_correlation = 1000  # need to be determined
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
                             'z',
                             'label_nuclear',
                             'centroid_nuclear',
                             'area_nuclear',
                             'mean_intensity_MYC_DNAFISH_in_nucleus',
                             'total_intensity_MYC_DNAFISH_in_nucleus'])

# LOAD IMAGE
for fov in range(total_fov):
    im_stack = skio.imread("%s%s_%s_%s.tif" % (master_folder, prefix, sample, fov + 1), plugin="tifffile")
    nuclear_seg_convex = np.zeros(shape=(im_stack.shape[0], im_stack.shape[2], im_stack.shape[3]))
    for z in range(im_stack.shape[0]):
        print("Start analyzing FOV %s/%s, z %s/%s" % (fov+1, total_fov, z+1, im_stack.shape[0]))
        img_nuclear = im_stack[z][1]
        img_FISH = im_stack[z][0]
        img_nuclear_seg = seg.nuclear_seg1(img_nuclear, local_factor=299, min_size=min_size_nuclear,
                                           max_size=max_size_nuclear)
        img_nuclear_seg_convex = seg.obj_to_convex(img_nuclear_seg)
        nuclear_seg_convex[z] = img_nuclear_seg_convex
        FISH_props = regionprops(img_nuclear_seg_convex, img_FISH)
        nuclear_props = regionprops(img_nuclear_seg_convex, img_nuclear)

        for i in range(len(nuclear_props)):
            # nuclear morphology
            nuclear_label = FISH_props[i].label
            nuclear_centroid = FISH_props[i].centroid
            nuclear_area = FISH_props[i].area

            # ecDNA related
            FISH_mean_intensity_nuclear = FISH_props[i].mean_intensity
            FISH_total_intensity_nuclear = nuclear_area * FISH_mean_intensity_nuclear

            data.loc[len(data.index)] = [fov, z, nuclear_label, nuclear_centroid, nuclear_area,
                                         FISH_mean_intensity_nuclear, FISH_total_intensity_nuclear]

    img.save_3d_image(nuclear_seg_convex, master_folder, 'nuclear_seg_convex_%s_%s.tiff' % (sample, i+1))

data.to_csv('%s%s_intensity.txt' % (master_folder, sample), index=False, sep='\t')


"""nuclear_props = regionprops(img_nuclear_seg)

        for i in range(len(nuclear_props)):
            print("Start image analysis z %s/%s, nuclear %s/%s..." % (
            z + 1, imstack.shape[0] + 1, i + 1, len(nuclear_props)))
            original_centroid = nuclear_props[i].centroid

            position = img.img_local_position(img_nuclear_seg, original_centroid, local_size)
            nuclear_seg = img.img_local_seg(img_nuclear_seg, position, i + 1)
            nuclear = img_nuclear[position[0]:position[1], position[2]:position[3]]
            centroid = regionprops(label(nuclear_seg))[0].centroid
            FISH = img_FISH[position[0]:position[1], position[2]:position[3]]

            nuclear_convex_local = regionprops(label(nuclear_seg))[0].convex_image
            centroid_convex = regionprops(label(nuclear_convex_local))[0].centroid
            nuclear_convex = img.image_paste(FISH, nuclear_convex_local, [int(centroid[0] - centroid_convex[0]),
                                                                          int(centroid[1] - centroid_convex[1])])

            nuclear_area = regionprops(label(nuclear_convex))[0].area
            FISH_mean_int = regionprops(label(nuclear_convex), FISH)[0].mean_intensity
            nuclear_circ = (4 * math.pi * nuclear_area) / (regionprops(label(nuclear_convex))[0].perimeter ** 2)

            _, r, g_FISH, dg_FISH = mat.auto_correlation(FISH, nuclear_seg, rmax)
            _, r, g_nuclear, dg_nuclear = mat.auto_correlation(nuclear, nuclear_seg, rmax)
            data.loc[len(data.index)] = \
                [3, z, i + 1, original_centroid, nuclear_area, nuclear_circ, FISH_mean_int, g_FISH, dg_FISH, g_nuclear,
                 dg_nuclear]

    data.to_csv('%sauto_correlation_%s.txt' % (master_folder, sample), index=False, sep='\t')














    for i in range(len(FISH_props)):
        # nuclear morphology
        nuclear_label = FISH_props[i].label
        nuclear_centroid = FISH_props[i].centroid
        nuclear_area = FISH_props[i].area

        # ecDNA related
        FISH_mean_intensity_nuclear = FISH_props[i].mean_intensity
        FISH_total_intensity_nuclear = nuclear_area * FISH_mean_intensity_nuclear

        # auto-correlation
        position = img.img_local_position(img_nuclear_seg_convex, nuclear_centroid, local_size)
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
        random_dot = random.choices(vector, cum_weights=vector_cum_weight, k=k_ddots)
        img_dot = np.zeros_like(nuclear_seg)
        for j in random_dot:
            img_dot[j[0]][j[1]] = img_dot[j[0]][j[1]] + 1

        _, r, g, dg = mat.auto_correlation(img_dot, nuclear_seg, rmax)

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
            [fov, nuclear_label, nuclear_centroid, nuclear_area, nuclear_major_axis, nuclear_minor_axis,
             nuclear_axis_ratio, nuclear_circularity, nuclear_eccentricity, nuclear_mean_intensity,
             nuclear_total_intensity, FISH_mean_intensity_nuclear,
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

data.to_csv('%s%s.txt' % (master_folder, sample), index=False, sep='\t')"""

print("DONE!")




