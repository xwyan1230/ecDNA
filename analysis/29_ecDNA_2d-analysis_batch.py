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
import shared.dataframe as dat

# input parameters
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220407_sp8_DMandHSR/HSR_singleZ/"
prefix = '20220407_DMandHSR_HSR_singleZ'
total_fov = 6
sample = 'HSR'
local_size = 150
rmax = 100
int_thresh_auto_correlation = 10000
radial_interval = 1
radial_max = 120
relative_radial_interval = 0.01
average_image_size = 400


data = pd.DataFrame(columns=['FOV', 'nuclear_label', 'nuclear_centroid', 'nuclear_area', 'nuclear_major_axis',
                             'nuclear_minor_axis', 'nuclear_axis_ratio', 'nuclear_circularity',
                             'nuclear_eccentricity', 'nuclear_mean_intensity', 'nuclear_total_intensity',
                             'FISH_mean_intensity_nuclear', 'FISH_total_intensity_nuclear',
                             'ecDNA_number', 'ecDNA_area', 'ecDNA_total_area', 'area_ratio', 'ecDNA_mean_area',
                             'ecDNA_max_area', 'ecDNA_mean_int',
                             'ecDNA_intensity', 'ecDNA_total_intensity', 'ecDNA_participating_coefficient',
                             'ecDNA_centroid', 'ecDNA_localization_from_centroid', 'ecDNA_distance_from_centroid',
                             'ecDNA_distance_from_edge',
                             'MYC_mean_intensity', 'MYC_total_intensity', 'g', 'dg', 'radial_distribution_from_centroid',
                             'radial_distribution_from_edge', 'radial_distribution_relative_r'])
FISH_sum = np.zeros(shape=(average_image_size, average_image_size))
nuclear_sum = np.zeros(shape=(average_image_size, average_image_size))
nuclear_seg_sum = np.zeros(shape=(average_image_size, average_image_size))
center = [average_image_size/2, average_image_size/2]

for fov in range(total_fov):
    print("Start analyzing FOV %s/%s" % (fov+1, total_fov))

    # load images
    img_nuclear = plt.imread('%s/%s_s%s_ch00.tif' % (master_folder, prefix, fov), format=None)
    img_FISH = plt.imread('%s/%s_s%s_ch02.tif' % (master_folder, prefix, fov), format=None)
    img_MYC = plt.imread('%s/%s_s%s_ch01.tif' % (master_folder, prefix, fov), format=None)

    # nuclear segmentation
    print("Start segmentation...")
    img_nuclear_seg = seg.nuclear_seg1(img_nuclear, local_factor=99, min_size=9000, max_size=50000)
    img_nuclear_seg = seg.filter_solidity(img_nuclear_seg)
    img_nuclear_seg_convex = seg.obj_to_convex(img_nuclear_seg)
    img_nuclear_seg_convex = dilation(img_nuclear_seg_convex, disk(3))

    # ecDNA segmentation
    _, img_FISH_seg = seg.find_organelle(img_FISH, 'na', extreme_val=18000, bg_val=seg.get_bg_int([img_FISH])[0],
                                         min_size=5, max_size=50000)
    img_FISH_seg = binary_dilation(img_FISH_seg, disk(10))
    for i in range(10):
        img_FISH_seg = binary_erosion(img_FISH_seg)
    img_FISH_seg = ndimage.binary_fill_holes(img_FISH_seg)
    img_FISH_seg = obj.remove_small(img_FISH_seg, min_size=20)

    # filter
    img_nuclear_seg_convex = seg.filter_mean_int(img_nuclear_seg_convex, img_FISH_seg, img_FISH, 20000)
    img_FISH_seg[img_nuclear_seg_convex == 0] = 0

    # props
    print("Start analyzing features...")

    FISH_props = regionprops(img_nuclear_seg_convex, img_FISH)
    nuclear_props = regionprops(img_nuclear_seg_convex, img_nuclear)
    MYC_props = regionprops(img_nuclear_seg_convex, img_MYC)

    for i in range(len(FISH_props)):
        # nuclear morphology
        nuclear_label = FISH_props[i].label
        nuclear_centroid = FISH_props[i].centroid
        nuclear_area = FISH_props[i].area
        nuclear_major_axis = FISH_props[i].major_axis_length
        nuclear_minor_axis = FISH_props[i].minor_axis_length
        nuclear_axis_ratio = nuclear_major_axis*1.0/nuclear_minor_axis
        nuclear_circularity = (4 * math.pi * nuclear_area) / (FISH_props[i].perimeter ** 2)
        nuclear_eccentricity = FISH_props[i].eccentricity
        nuclear_mean_intensity = nuclear_props[i].mean_intensity
        nuclear_total_intensity = nuclear_area * nuclear_mean_intensity

        # ecDNA related
        FISH_mean_intensity_nuclear = FISH_props[i].mean_intensity
        FISH_total_intensity_nuclear = nuclear_area * FISH_mean_intensity_nuclear
        img_FISH_seg_temp = img_FISH_seg.copy()
        img_FISH_seg_temp[img_nuclear_seg_convex != FISH_props[i].label] = 0
        ecDNA_props = regionprops(label(img_FISH_seg_temp), img_FISH)
        ecDNA_number = len(ecDNA_props)
        ecDNA_area = [ecDNA_props[j].area for j in range(ecDNA_number)]
        ecDNA_total_area = sum(ecDNA_area)
        area_ratio = ecDNA_total_area/nuclear_area
        ecDNA_mean_area = ecDNA_total_area/ecDNA_number
        ecDNA_max_area = max(ecDNA_area)
        ecDNA_mean_int = [ecDNA_props[j].mean_intensity for j in range(ecDNA_number)]
        ecDNA_intensity = [ecDNA_mean_int[j]*ecDNA_area[j] for j in range(ecDNA_number)]
        ecDNA_total_intensity = sum(ecDNA_intensity)
        ecDNA_participating_coefficient = ecDNA_total_intensity * 1.0 / FISH_total_intensity_nuclear
        ecDNA_centroid = [ecDNA_props[j].centroid for j in range(ecDNA_number)]
        ecDNA_localization_from_centroid = [tuple(np.array(ecDNA_centroid[j])-np.array(nuclear_centroid))
                                            for j in range(ecDNA_number)]
        ecDNA_distance_from_centroid = \
            [(ecDNA_localization_from_centroid[j][0]**2+ecDNA_localization_from_centroid[j][1]**2)**0.5
             for j in range(ecDNA_number)]
        _, distance_map = medial_axis(img_nuclear_seg_convex, return_distance=True)
        ecDNA_distance_from_edge = [distance_map[round(ecDNA_centroid[j][0])][round(ecDNA_centroid[j][1])]
                                    for j in range(ecDNA_number)]

        # MYC expression
        MYC_mean_intensity = MYC_props[i].mean_intensity
        MYC_total_intensity = nuclear_area * MYC_mean_intensity

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

        """
        # old algorithm for radial distribution analysis
        for m in range(len(nuclear_seg)):
            for n in range(len(nuclear_seg[0])):
                if nuclear_seg[m][n] != 0:
                    pixel_localization_from_centroid = tuple(np.array([m, n]) - np.array(local_nuclear_centroid))
                    pixel_distance_from_centroid = (pixel_localization_from_centroid[0]**2 +
                                                    pixel_localization_from_centroid[1]**2)**0.5
                    pixel_distance_from_edge = local_distance_map[m][n]
                    pixel_relative_r = pixel_distance_from_centroid/(pixel_distance_from_edge +
                                                                     pixel_distance_from_centroid)
                    pixel_FISH_intensity = FISH[m][n]
                    pixel.loc[len(pixel.index)] = [pixel_distance_from_centroid, pixel_distance_from_edge,
                                                   pixel_relative_r, pixel_FISH_intensity]

        radial_distribution_from_centroid = \
            dat.radial_distribution(pixel, 'pixel_distance_from_centroid', 'pixel_FISH_intensity', radial_interval,
                                    radial_max)
        radial_distribution_from_edge = \
            dat.radial_distribution(pixel, 'pixel_distance_from_edge', 'pixel_FISH_intensity', radial_interval,
                                    radial_max)
        radial_distribution_relative_r = \
            dat.radial_distribution(pixel, 'pixel_relative_r', 'pixel_FISH_intensity', relative_radial_interval, 1)"""

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
             g, dg, radial_distribution_from_centroid, radial_distribution_from_edge, radial_distribution_relative_r]

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




