import matplotlib.pyplot as plt
from skimage.measure import label, regionprops
import shared.segmentation as seg
from skimage.morphology import binary_dilation, binary_erosion, disk, dilation
import shared.objects as obj
from scipy import ndimage
import pandas as pd
import math
import numpy as np

# input parameters
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220407_sp8_DMandHSR/DM_singleZ/"
prefix = '20220407_DMandHSR_DM_singleZ'
total_fov = 6
sample = 'DM'

for fov in range(total_fov):
    print("Start analyzing FOV %s/%s" % (fov+1, total_fov))

    # load images
    img_nuclear = plt.imread('%s/%s_s%s_ch00.tif' % (master_folder, prefix, fov), format=None)
    img_FISH = plt.imread('%s/%s_s%s_ch02.tif' % (master_folder, prefix, fov), format=None)

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
                                 'ecDNA_mean_area', 'ecDNA_max_area', 'ecDNA_mean_int', 'ecDNA_intensity',
                                 'ecDNA_total_intensity', 'ecDNA_participating_coefficient', 'ecDNA_centroid',
                                 'ecDNA_localization_from_centroid', 'ecDNA_distance_from_centroid'])

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
        ecDNA_mean_area = sum(ecDNA_area)/ecDNA_number
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

        data.loc[len(data.index)] = \
            [fov, nuclear_label, nuclear_centroid, nuclear_area, nuclear_major_axis, nuclear_minor_axis, nuclear_axis_ratio,
             nuclear_circularity, nuclear_eccentricity, nuclear_FISH_mean_intensity, nuclear_total_intensity, ecDNA_number,
             ecDNA_area, ecDNA_mean_area, ecDNA_max_area, ecDNA_mean_int, ecDNA_intensity, ecDNA_total_intensity,
             ecDNA_participating_coefficient, ecDNA_centroid, ecDNA_localization_from_centroid, ecDNA_distance_from_centroid]

data.to_csv('%s%s.txt' % (master_folder, sample), index=False, sep='\t')

print("DONE!")




