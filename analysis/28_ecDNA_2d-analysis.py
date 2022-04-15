import matplotlib.pyplot as plt
from skimage.measure import label, regionprops
import shared.segmentation as seg
from skimage.morphology import binary_dilation, binary_erosion, disk, dilation
import napari
import shared.objects as obj
from scipy import ndimage
import pandas as pd
import math
import numpy as np
import shared.image as img
import random
import shared.math as mat

# input parameters
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220407_sp8_DMandHSR/DM_singleZ/"
prefix = '20220407_DMandHSR_DM_singleZ'
fov = 0
local_size = 150
dead_int = 10000
rmax = 100

# load images
img_nuclear = plt.imread('%s/%s_s0_ch00.tif' % (master_folder, prefix), format=None)
img_FISH = plt.imread('%s/%s_s0_ch02.tif' % (master_folder, prefix), format=None)
img_MYC = plt.imread('%s/%s_s0_ch01.tif' % (master_folder, prefix), format=None)

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
data = pd.DataFrame(columns=['FOV', 'nuclear_label', 'nuclear_centroid', 'nuclear_area', 'nuclear_major_axis',
                             'nuclear_minor_axis', 'nuclear_axis_ratio', 'nuclear_circularity', 'nuclear_eccentricity',
                             'nuclear_FISH_mean_intensity', 'nuclear_total_intensity', 'ecDNA_number', 'ecDNA_area',
                             'ecDNA_total_area', 'area_ratio',
                             'ecDNA_mean_area', 'ecDNA_max_area', 'ecDNA_mean_int', 'ecDNA_intensity',
                             'ecDNA_total_intensity', 'ecDNA_participating_coefficient', 'ecDNA_centroid',
                             'ecDNA_localization_from_centroid', 'ecDNA_distance_from_centroid', 'MYC_mean_intensity',
                             'MYC_total_intensity', 'g', 'dg'])

props = regionprops(img_nuclear_seg_convex, img_FISH)
for i in range(len(props)):
    nuclear_label = props[i].label
    nuclear_centroid = props[i].centroid
    nuclear_area = props[i].area
    nuclear_major_axis = props[i].major_axis_length
    nuclear_minor_axis = props[i].minor_axis_length
    nuclear_axis_ratio = nuclear_major_axis*1.0/nuclear_minor_axis
    nuclear_circularity = (4 * math.pi * nuclear_area) / (props[i].perimeter ** 2)
    nuclear_eccentricity = props[i].eccentricity
    nuclear_FISH_mean_intensity = props[i].mean_intensity
    nuclear_total_intensity = nuclear_area * nuclear_FISH_mean_intensity

    img_FISH_seg_temp = img_FISH_seg.copy()
    img_FISH_seg_temp[img_nuclear_seg_convex != props[i].label] = 0
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
    ecDNA_participating_coefficient = ecDNA_total_intensity*1.0/nuclear_total_intensity
    ecDNA_centroid = [ecDNA_props[j].centroid for j in range(ecDNA_number)]
    ecDNA_localization_from_centroid = [tuple(np.array(ecDNA_centroid[j])-np.array(nuclear_centroid))
                                        for j in range(ecDNA_number)]
    ecDNA_distance_from_centroid = \
        [(ecDNA_localization_from_centroid[j][0]**2+ecDNA_localization_from_centroid[j][1]**2)**0.5
         for j in range(ecDNA_number)]

    MYC_props = regionprops(img_nuclear_seg_convex, img_MYC)
    MYC_mean_intensity = MYC_props[i].mean_intensity
    MYC_total_intensity = nuclear_area*MYC_mean_intensity

    position = img.img_local_position(img_nuclear_seg_convex, nuclear_centroid, local_size)
    nuclear_seg = img.img_local_seg(img_nuclear_seg_convex, position, i + 1)
    FISH = img_FISH[position[0]:position[1], position[2]:position[3]]
    vector = []
    vector_cum_weight = []
    weight = 0
    for m in range(len(nuclear_seg)):
        for n in range(len(nuclear_seg[0])):
            if nuclear_seg[m][n] == 1:
                vector.append([m, n])
                if FISH[m][n] > dead_int:
                    weight = weight + FISH[m][n] - dead_int
                vector_cum_weight.append(weight)
    random_dot = random.choices(vector, cum_weights=vector_cum_weight, k=100)
    img_dot = np.zeros_like(nuclear_seg)
    for j in random_dot:
        img_dot[j[0]][j[1]] = img_dot[j[0]][j[1]] + 1

    _, r, g, dg = mat.auto_correlation(img_dot, nuclear_seg, rmax)

    """viewer = napari.Viewer()
    viewer.add_image(nuclear_seg)
    viewer.add_image(FISH)
    napari.run()
    viewer = napari.Viewer()
    viewer.add_image(img_dot)
    napari.run()"""

    data.loc[len(data.index)] = \
        [fov, nuclear_label, nuclear_centroid, nuclear_area, nuclear_major_axis, nuclear_minor_axis, nuclear_axis_ratio,
         nuclear_circularity, nuclear_eccentricity, nuclear_FISH_mean_intensity, nuclear_total_intensity, ecDNA_number,
         ecDNA_area, ecDNA_total_area, area_ratio, ecDNA_mean_area, ecDNA_max_area, ecDNA_mean_int, ecDNA_intensity,
         ecDNA_total_intensity,
         ecDNA_participating_coefficient, ecDNA_centroid, ecDNA_localization_from_centroid, ecDNA_distance_from_centroid,
         MYC_mean_intensity, MYC_total_intensity, g, dg]

#data.to_csv('%stest.txt' % master_folder, index=False, sep='\t')

viewer = napari.Viewer()
viewer.add_image(img_nuclear_seg_convex, colormap='viridis', blending='additive')
viewer.add_image(img_FISH_seg, blending='additive')
napari.run()



